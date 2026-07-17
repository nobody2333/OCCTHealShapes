#include "HealClosedNarrowGap.h"

#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>

namespace {
constexpr int kMaxAutomaticRepairIterations = 300;

struct BoundaryEdgeUse {
    TopoDS_Edge edge;
    TopoDS_Face ownerFace;
    TopoDS_Vertex firstVertex;
    TopoDS_Vertex lastVertex;
};

struct BoundaryLoop {
    std::vector<TopoDS_Edge> edges;
    std::vector<TopoDS_Face> ownerFaces;
    std::vector<TopoDS_Vertex> vertices;
};

double clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

std::string failureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
}

int listExtent(const TopTools_ListOfShape& shapes)
{
    int count = 0;
    for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
        ++count;
    }
    return count;
}

std::string describeEdgeFaceCounts(
    const TopoDS_Shape& shape,
    const std::vector<TopoDS_Edge>& edges)
{
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);

    std::ostringstream os;
    os << "edgeCount=" << edges.size() << ", adjacentFaceCounts=[";
    for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
        if (edgeIndex > 0) {
            os << ", ";
        }

        int faceCount = 0;
        if (!edges[edgeIndex].IsNull() && edgeFaceMap.Contains(edges[edgeIndex])) {
            faceCount = listExtent(edgeFaceMap.FindFromKey(edges[edgeIndex]));
        }
        os << edgeIndex << ":" << faceCount;
    }
    os << "]";
    return os.str();
}

TopoDS_Face firstFaceInList(const TopTools_ListOfShape& shapes)
{
    for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (!face.IsNull()) {
            return face;
        }
    }
    return TopoDS_Face();
}

bool getEdgeVertices(const TopoDS_Edge& edge, TopoDS_Vertex& first, TopoDS_Vertex& last)
{
    first.Nullify();
    last.Nullify();
    if (edge.IsNull()) {
        return false;
    }

    TopExp::Vertices(edge, first, last, Standard_True);
    return !first.IsNull() && !last.IsNull();
}

TopoDS_Edge orientEdgeFromTo(const TopoDS_Edge& edge, const TopoDS_Vertex& from, const TopoDS_Vertex& to)
{
    TopoDS_Edge oriented = edge;
    if (edge.IsNull() || from.IsNull() || to.IsNull()) {
        return oriented;
    }

    TopoDS_Vertex first;
    TopoDS_Vertex last;
    if (!getEdgeVertices(oriented, first, last)) {
        return oriented;
    }

    if (first.IsSame(from) && last.IsSame(to)) {
        return oriented;
    }

    if (first.IsSame(to) && last.IsSame(from)) {
        oriented.Reverse();
    }
    return oriented;
}

bool edgeHasVertex(const BoundaryEdgeUse& edgeUse, const TopoDS_Vertex& vertex)
{
    return !vertex.IsNull()
        && (edgeUse.firstVertex.IsSame(vertex) || edgeUse.lastVertex.IsSame(vertex));
}

TopoDS_Vertex otherVertexOf(const BoundaryEdgeUse& edgeUse, const TopoDS_Vertex& vertex)
{
    if (edgeUse.firstVertex.IsSame(vertex)) {
        return edgeUse.lastVertex;
    }
    if (edgeUse.lastVertex.IsSame(vertex)) {
        return edgeUse.firstVertex;
    }
    return TopoDS_Vertex();
}

bool isSameFaceInVector(const std::vector<TopoDS_Face>& faces, const TopoDS_Face& face)
{
    for (const TopoDS_Face& existingFace : faces) {
        if (existingFace.IsSame(face)) {
            return true;
        }
    }
    return false;
}

bool isBoundaryEdgeUsable(const TopoDS_Edge& edge, const TopoDS_Face& ownerFace)
{
    if (edge.IsNull() || ownerFace.IsNull()) {
        return false;
    }

    try {
        if (BRep_Tool::Degenerated(edge) == Standard_True) {
            return false;
        }
        if (BRep_Tool::IsClosed(edge) || BRep_Tool::IsClosed(edge, ownerFace)) {
            return false;
        }
    }
    catch (Standard_Failure&) {
        return false;
    }

    return GeometryUtils::ComputeEdgeLength(edge) > Precision::Confusion();
}

std::vector<BoundaryLoop> buildBoundaryLoops(const TopoDS_Shape& shape)
{
    std::vector<BoundaryLoop> loops;
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);

    std::vector<BoundaryEdgeUse> boundaryEdges;
    for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
        if (listExtent(edgeFaceMap.FindFromIndex(i)) != 1) {
            continue;
        }

        const TopoDS_Edge edge = TopoDS::Edge(edgeFaceMap.FindKey(i));
        const TopoDS_Face ownerFace = firstFaceInList(edgeFaceMap.FindFromIndex(i));
        if (!isBoundaryEdgeUsable(edge, ownerFace)) {
            continue;
        }

        BoundaryEdgeUse edgeUse;
        edgeUse.edge = edge;
        edgeUse.ownerFace = ownerFace;
        if (!getEdgeVertices(edge, edgeUse.firstVertex, edgeUse.lastVertex)) {
            continue;
        }
        boundaryEdges.push_back(edgeUse);
    }

    std::vector<bool> consumed(boundaryEdges.size(), false);
    for (std::size_t startIndex = 0; startIndex < boundaryEdges.size(); ++startIndex) {
        if (consumed[startIndex]) {
            continue;
        }

        BoundaryLoop loop;
        const TopoDS_Vertex loopStart = boundaryEdges[startIndex].firstVertex;
        TopoDS_Vertex currentVertex = boundaryEdges[startIndex].lastVertex;
        loop.vertices.push_back(loopStart);
        loop.vertices.push_back(currentVertex);
        loop.edges.push_back(orientEdgeFromTo(boundaryEdges[startIndex].edge, loopStart, currentVertex));
        loop.ownerFaces.push_back(boundaryEdges[startIndex].ownerFace);
        consumed[startIndex] = true;

        bool closed = false;
        bool failed = false;
        while (!failed && loop.edges.size() <= boundaryEdges.size()) {
            if (currentVertex.IsSame(loopStart)) {
                closed = true;
                break;
            }

            int nextIndex = -1;
            for (std::size_t i = 0; i < boundaryEdges.size(); ++i) {
                if (consumed[i] || !edgeHasVertex(boundaryEdges[i], currentVertex)) {
                    continue;
                }
                if (nextIndex >= 0) {
                    failed = true;
                    break;
                }
                nextIndex = static_cast<int>(i);
            }

            if (failed || nextIndex < 0) {
                failed = true;
                break;
            }

            const BoundaryEdgeUse& nextUse = boundaryEdges[static_cast<std::size_t>(nextIndex)];
            TopoDS_Vertex nextVertex = otherVertexOf(nextUse, currentVertex);
            if (nextVertex.IsNull()) {
                failed = true;
                break;
            }

            loop.edges.push_back(orientEdgeFromTo(nextUse.edge, currentVertex, nextVertex));
            loop.ownerFaces.push_back(nextUse.ownerFace);
            loop.vertices.push_back(nextVertex);
            consumed[static_cast<std::size_t>(nextIndex)] = true;
            currentVertex = nextVertex;
        }

        if (closed && loop.edges.size() >= 4 && loop.vertices.size() == loop.edges.size() + 1) {
            loops.push_back(loop);
        }
    }

    return loops;
}
} // namespace

HealClosedNarrowGap::HealClosedNarrowGap(const Params& params)
    : params_(params)
{
    params_.maxGapWidth = std::max(params_.maxGapWidth, Precision::Confusion());
    params_.minLengthRatio = clamp01(params_.minLengthRatio);
    params_.minLengthWidthRatio = std::max(1.0, params_.minLengthWidthRatio);
    params_.sampleCount = std::max(2, params_.sampleCount);
    params_.minSampleCount = std::max(2, params_.minSampleCount);
    params_.maxSampleCount = std::max(params_.minSampleCount, params_.maxSampleCount);
    if (params_.targetSampleSpacing <= Precision::Confusion()) {
        params_.useAdaptiveSampleCount = false;
    }
    params_.midlineRepairTolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
}

void HealClosedNarrowGap::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = shape;
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
}

int HealClosedNarrowGap::perform()
{
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    TopoDS_Shape currentShape = inputShape_;
    const int iterationLimit = computeAutomaticIterationLimit(currentShape);
    {
        std::ostringstream os;
        os << "Automatic repair iteration limit = " << iterationLimit;
        appendMessage(os.str());
    }
    {
        std::ostringstream os;
        os << "Parameters: maxGapWidth=" << params_.maxGapWidth
           << ", minLengthRatio=" << params_.minLengthRatio
           << ", minLengthWidthRatio=" << params_.minLengthWidthRatio
           << ", maxTipAngle=" << params_.maxTipAngle
           << ", sampleCount=" << params_.sampleCount
           << ", useAdaptiveSampleCount="
           << (params_.useAdaptiveSampleCount ? "true" : "false")
           << ", minSampleCount=" << params_.minSampleCount
           << ", maxSampleCount=" << params_.maxSampleCount
           << ", targetSampleSpacing=" << params_.targetSampleSpacing
           << ", useMidlineCollapse=" << (params_.useMidlineCollapse ? "true" : "false")
           << ", fallbackToSideReplacement=" << (params_.fallbackToSideReplacement ? "true" : "false")
           << ", midlineRepairTolerance=" << params_.midlineRepairTolerance
           << ", runShapeFixAfterRepair="
           << (params_.runShapeFixAfterRepair ? "true" : "false");
        appendMessage(os.str());
    }

    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Original shape", currentShape);

    for (int iter = 0; iter < iterationLimit; ++iter) {
        RepairDiagnostics::AppendIterationHeader(messages_, iter + 1, iter > 0);

        std::vector<ClosedGapCandidate> candidates = extractCandidates(currentShape);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const ClosedGapCandidate& lhs, const ClosedGapCandidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "Candidate count = " << candidates.size();
            appendMessage(os.str());
        }

        if (candidates.empty()) {
            break;
        }

        bool repairedThisIteration = false;
        for (const ClosedGapCandidate& candidate : candidates) {
            TopoDS_Shape repairedShape;
            std::string errorMessage;
            if (!tryRepairCandidate(currentShape, candidate, repairedShape, errorMessage)) {
                std::ostringstream os;
                os << "Candidate failed:"
                   << "\n  chainAEdges=" << candidate.chainA.edges.size()
                   << "\n  chainBEdges=" << candidate.chainB.edges.size()
                   << "\n  maxWidth=" << candidate.maxWidth
                   << "\n  meanWidth=" << candidate.meanWidth
                   << "\n  lengthRatio=" << candidate.lengthRatio
                   << "\n  lengthWidthRatio=" << candidate.lengthWidthRatio
                   << "\n  endpointSharpness=" << candidate.endpointSharpness
                   << "\n  confidence=" << candidate.confidence
                   << "\n  reason:"
                   << "\n" << RepairDiagnostics::IndentMultiline(
                       RepairDiagnostics::MakeFailureReasonReadable(errorMessage),
                       "    ");
                appendMessage(os.str());
                continue;
            }

            currentShape = repairedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++repairedCount_;
            repairedThisIteration = true;

            std::ostringstream os;
            os << "Repaired candidate: chainAEdges=" << candidate.chainA.edges.size()
               << ", chainBEdges=" << candidate.chainB.edges.size()
               << ", maxWidth=" << candidate.maxWidth
               << ", meanWidth=" << candidate.meanWidth
               << ", lengthRatio=" << candidate.lengthRatio
               << ", lengthWidthRatio=" << candidate.lengthWidthRatio
               << ", endpointSharpness=" << candidate.endpointSharpness
               << ", confidence=" << candidate.confidence;
            appendMessage(os.str());
            break;
        }

        if (!repairedThisIteration) {
            appendMessage("Repair stopped: no candidate was submitted in this iteration");
            break;
        }

        if (iter + 1 == iterationLimit) {
            appendMessage("Repair stopped: automatic iteration limit reached");
        }
    }

    resultShape_ = currentShape;
    appendMessage("");
    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Result shape", resultShape_);
    return repairedCount_;
}

int HealClosedNarrowGap::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

std::vector<HealClosedNarrowGap::ClosedGapCandidate> HealClosedNarrowGap::extractCandidates(
    const TopoDS_Shape& shape) const
{
    std::vector<ClosedGapCandidate> candidates;
    if (shape.IsNull()) {
        return candidates;
    }

    const std::vector<BoundaryLoop> loops = buildBoundaryLoops(shape);
    for (const BoundaryLoop& loop : loops) {
        const int edgeCount = static_cast<int>(loop.edges.size());
        if (edgeCount < 4) {
            continue;
        }

        std::vector<int> tipIndices;
        std::vector<double> tipScores(edgeCount, 0.0);
        for (int i = 0; i < edgeCount; ++i) {
            const int prevIndex = (i + edgeCount - 1) % edgeCount;
            const double angle = GeometryUtils::ComputeIntersectionAngle(
                loop.edges[prevIndex],
                loop.edges[i],
                loop.vertices[static_cast<std::size_t>(i)]);
            if (angle < 0.0 || angle > params_.maxTipAngle) {
                continue;
            }

            tipIndices.push_back(i);
            tipScores[static_cast<std::size_t>(i)] =
                1.0 - clamp01(angle / std::max(params_.maxTipAngle, Precision::Confusion()));
        }

        if (tipIndices.size() < 2) {
            continue;
        }

        for (std::size_t lhs = 0; lhs < tipIndices.size(); ++lhs) {
            for (std::size_t rhs = lhs + 1; rhs < tipIndices.size(); ++rhs) {
                const int startIndex = tipIndices[lhs];
                const int endIndex = tipIndices[rhs];

                EdgeChain chainA;
                EdgeChain pathB;
                int guard = 0;
                int index = startIndex;
                chainA.vertices.push_back(loop.vertices[static_cast<std::size_t>(startIndex)]);
                while (index != endIndex && guard++ <= edgeCount) {
                    const int nextIndex = (index + 1) % edgeCount;
                    const TopoDS_Vertex from = loop.vertices[static_cast<std::size_t>(index)];
                    const TopoDS_Vertex to = loop.vertices[static_cast<std::size_t>(nextIndex)];
                    const TopoDS_Edge edge = orientEdgeFromTo(loop.edges[static_cast<std::size_t>(index)], from, to);
                    chainA.edges.push_back(edge);
                    chainA.ownerFaces.push_back(loop.ownerFaces[static_cast<std::size_t>(index)]);
                    chainA.vertices.push_back(to);
                    chainA.length += GeometryUtils::ComputeEdgeLength(edge);
                    index = nextIndex;
                }

                guard = 0;
                index = endIndex;
                pathB.vertices.push_back(loop.vertices[static_cast<std::size_t>(endIndex)]);
                while (index != startIndex && guard++ <= edgeCount) {
                    const int nextIndex = (index + 1) % edgeCount;
                    const TopoDS_Vertex from = loop.vertices[static_cast<std::size_t>(index)];
                    const TopoDS_Vertex to = loop.vertices[static_cast<std::size_t>(nextIndex)];
                    const TopoDS_Edge edge = orientEdgeFromTo(loop.edges[static_cast<std::size_t>(index)], from, to);
                    pathB.edges.push_back(edge);
                    pathB.ownerFaces.push_back(loop.ownerFaces[static_cast<std::size_t>(index)]);
                    pathB.vertices.push_back(to);
                    pathB.length += GeometryUtils::ComputeEdgeLength(edge);
                    index = nextIndex;
                }

                if (chainA.edges.empty() || pathB.edges.empty()) {
                    continue;
                }

                EdgeChain chainB;
                chainB.vertices.push_back(pathB.vertices.back());
                for (int i = static_cast<int>(pathB.edges.size()) - 1; i >= 0; --i) {
                    const TopoDS_Vertex from = pathB.vertices[static_cast<std::size_t>(i + 1)];
                    const TopoDS_Vertex to = pathB.vertices[static_cast<std::size_t>(i)];
                    const TopoDS_Edge edge = orientEdgeFromTo(pathB.edges[static_cast<std::size_t>(i)], from, to);
                    chainB.edges.push_back(edge);
                    chainB.ownerFaces.push_back(pathB.ownerFaces[static_cast<std::size_t>(i)]);
                    chainB.vertices.push_back(to);
                    chainB.length += GeometryUtils::ComputeEdgeLength(edge);
                }

                if (chainA.length <= Precision::Confusion() || chainB.length <= Precision::Confusion()) {
                    continue;
                }

                ClosedGapCandidate candidate;
                candidate.chainA = chainA;
                candidate.chainB = chainB;
                candidate.loopEdges = loop.edges;
                candidate.tipVertices.push_back(loop.vertices[static_cast<std::size_t>(startIndex)]);
                candidate.tipVertices.push_back(loop.vertices[static_cast<std::size_t>(endIndex)]);
                candidate.endpointSharpness =
                    0.5 * (tipScores[static_cast<std::size_t>(startIndex)]
                         + tipScores[static_cast<std::size_t>(endIndex)]);

                double gapSum = 0.0;
                bool samplingFailed = false;
                for (int sample = 0; sample < params_.sampleCount; ++sample) {
                    const double t = params_.sampleCount == 1
                        ? 0.0
                        : static_cast<double>(sample) / static_cast<double>(params_.sampleCount - 1);

                    gp_Pnt pointA;
                    gp_Pnt pointB;
                    auto sampleChainPoint = [](const EdgeChain& chain, double normalized, gp_Pnt& point) {
                        if (chain.edges.empty() || chain.vertices.size() != chain.edges.size() + 1 ||
                            chain.length <= Precision::Confusion()) {
                            return false;
                        }

                        const double targetLength = chain.length * clamp01(normalized);
                        double accumulated = 0.0;
                        for (std::size_t i = 0; i < chain.edges.size(); ++i) {
                            const TopoDS_Edge edge = chain.edges[i];
                            const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
                            if (edgeLength <= Precision::Confusion()) {
                                continue;
                            }

                            const bool isLastEdge = i + 1 == chain.edges.size();
                            if (!isLastEdge && accumulated + edgeLength < targetLength) {
                                accumulated += edgeLength;
                                continue;
                            }

                            const double local = isLastEdge
                                ? clamp01((targetLength - accumulated) / edgeLength)
                                : clamp01((targetLength - accumulated) / edgeLength);
                            const double firstParam = BRep_Tool::Parameter(chain.vertices[i], edge);
                            const double lastParam = BRep_Tool::Parameter(chain.vertices[i + 1], edge);
                            const double param = firstParam + (lastParam - firstParam) * local;

                            BRepAdaptor_Curve curve(edge);
                            point = curve.Value(param);
                            return true;
                        }

                        return false;
                    };

                    if (!sampleChainPoint(candidate.chainA, t, pointA) ||
                        !sampleChainPoint(candidate.chainB, t, pointB)) {
                        samplingFailed = true;
                        break;
                    }

                    const double width = pointA.Distance(pointB);
                    gapSum += width;
                    candidate.maxWidth = std::max(candidate.maxWidth, width);
                }

                if (samplingFailed) {
                    continue;
                }

                candidate.meanWidth = gapSum / static_cast<double>(params_.sampleCount);
                candidate.lengthRatio = std::min(chainA.length, chainB.length) / std::max(chainA.length, chainB.length);
                candidate.lengthWidthRatio =
                    0.5 * (chainA.length + chainB.length) / std::max(candidate.meanWidth, Precision::Confusion());

                if (candidate.maxWidth > params_.maxGapWidth ||
                    candidate.lengthRatio + Precision::Confusion() < params_.minLengthRatio ||
                    candidate.lengthWidthRatio + Precision::Confusion() < params_.minLengthWidthRatio) {
                    continue;
                }

                const double widthScore = 1.0 - clamp01(candidate.maxWidth / params_.maxGapWidth);
                const double ratioScore = clamp01(candidate.lengthRatio);
                const double slenderScore =
                    clamp01(candidate.lengthWidthRatio / std::max(params_.minLengthWidthRatio * 2.0, 1.0));
                candidate.confidence =
                    0.40 * widthScore +
                    0.25 * slenderScore +
                    0.20 * ratioScore +
                    0.15 * clamp01(candidate.endpointSharpness);
                candidates.push_back(candidate);
            }
        }
    }

    return candidates;
}

bool HealClosedNarrowGap::sampleChainPoint(const EdgeChain& chain, double normalized, gp_Pnt& point) const
{
    if (chain.edges.empty() || chain.vertices.size() != chain.edges.size() + 1 ||
        chain.length <= Precision::Confusion()) {
        return false;
    }

    const double targetLength = chain.length * clamp01(normalized);
    double accumulated = 0.0;
    for (std::size_t i = 0; i < chain.edges.size(); ++i) {
        const TopoDS_Edge edge = chain.edges[i];
        const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
        if (edgeLength <= Precision::Confusion()) {
            continue;
        }

        const bool isLastEdge = i + 1 == chain.edges.size();
        if (!isLastEdge && accumulated + edgeLength < targetLength) {
            accumulated += edgeLength;
            continue;
        }

        const double local = clamp01((targetLength - accumulated) / edgeLength);
        const double firstParam = BRep_Tool::Parameter(chain.vertices[i], edge);
        const double lastParam = BRep_Tool::Parameter(chain.vertices[i + 1], edge);
        const double param = firstParam + (lastParam - firstParam) * local;

        BRepAdaptor_Curve curve(edge);
        point = curve.Value(param);
        return true;
    }

    return false;
}

bool HealClosedNarrowGap::buildMidlineChain(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    if (candidate.chainA.edges.empty() || candidate.chainB.edges.empty()) {
        errorMessage = "candidate chains are empty";
        return false;
    }
    if (candidate.chainA.edges.size() != candidate.chainB.edges.size()) {
        errorMessage = "midline collapse requires one-to-one source edge correspondence";
        return false;
    }
    if (candidate.chainA.vertices.size() != candidate.chainA.edges.size() + 1 ||
        candidate.chainB.vertices.size() != candidate.chainB.edges.size() + 1) {
        errorMessage = "candidate chain vertex count is invalid";
        return false;
    }
    if (!candidate.chainA.vertices.front().IsSame(candidate.chainB.vertices.front()) ||
        !candidate.chainA.vertices.back().IsSame(candidate.chainB.vertices.back())) {
        errorMessage = "candidate chains do not share the same endpoints";
        return false;
    }

    try {
        const std::size_t edgeCount = candidate.chainA.edges.size();
        std::vector<TopoDS_Vertex> midVertices;
        midVertices.reserve(edgeCount + 1);

        BRep_Builder builder;
        for (std::size_t vertexIndex = 0; vertexIndex <= edgeCount; ++vertexIndex) {
            if (vertexIndex == 0) {
                // 两侧起始尖点已经拓扑重合，直接沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.front());
            }
            else if (vertexIndex == edgeCount) {
                // 两侧结束尖点同样沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.back());
            }
            else {
                const TopoDS_Vertex& vertexA = candidate.chainA.vertices[vertexIndex];
                const TopoDS_Vertex& vertexB = candidate.chainB.vertices[vertexIndex];
                if (vertexA.IsNull() || vertexB.IsNull()) {
                    errorMessage = "candidate contains a null corresponding internal vertex";
                    return false;
                }

                const gp_Pnt pointA = BRep_Tool::Pnt(vertexA);
                const gp_Pnt pointB = BRep_Tool::Pnt(vertexB);
                const double vertexGap = pointA.Distance(pointB);
                if (vertexGap > params_.maxGapWidth + params_.midlineRepairTolerance) {
                    errorMessage = "corresponding internal chain vertices are too far apart";
                    return false;
                }

                const gp_Pnt midpoint(
                    0.5 * (pointA.X() + pointB.X()),
                    0.5 * (pointA.Y() + pointB.Y()),
                    0.5 * (pointA.Z() + pointB.Z()));

                // 新顶点还会替换两侧旧顶点，因此容差必须覆盖旧顶点到中点的移动距离。
                const double vertexTolerance = std::max({
                    params_.midlineRepairTolerance,
                    BRep_Tool::Tolerance(vertexA),
                    BRep_Tool::Tolerance(vertexB),
                    midpoint.Distance(pointA) + Precision::Confusion(),
                    midpoint.Distance(pointB) + Precision::Confusion()
                });

                TopoDS_Vertex midVertex;
                builder.MakeVertex(midVertex, midpoint, vertexTolerance);
                if (midVertex.IsNull()) {
                    errorMessage = "failed to create midline vertex";
                    return false;
                }
                midVertices.push_back(midVertex);
            }
        }

        midlineData.chain.vertices = midVertices;
        midlineData.replacementsForChainA.resize(edgeCount);
        midlineData.replacementsForChainB.resize(edgeCount);

        auto sampleOrientedEdge = [&](const EdgeChain& chain,
                                      std::size_t edgeIndex,
                                      int sampleCount,
                                      std::vector<gp_Pnt>& sampledPoints,
                                      const char* chainName) {
            sampledPoints.clear();
            sampledPoints.reserve(static_cast<std::size_t>(sampleCount));

            const TopoDS_Edge& edge = chain.edges[edgeIndex];
            const TopoDS_Vertex& startVertex = chain.vertices[edgeIndex];
            const TopoDS_Vertex& endVertex = chain.vertices[edgeIndex + 1];
            const double startParameter = BRep_Tool::Parameter(startVertex, edge);
            const double endParameter = BRep_Tool::Parameter(endVertex, edge);
            if (!std::isfinite(startParameter) || !std::isfinite(endParameter) ||
                std::abs(endParameter - startParameter) <= Precision::PConfusion()) {
                errorMessage = std::string(chainName) + " edge has an invalid parameter range";
                return false;
            }

            BRepAdaptor_Curve curve(edge);
            const double lowerParameter = std::min(startParameter, endParameter);
            const double upperParameter = std::max(startParameter, endParameter);
            GCPnts_UniformAbscissa uniformSampler(
                curve,
                sampleCount,
                lowerParameter,
                upperParameter,
                Precision::Confusion());
            const bool useUniformAbscissa =
                uniformSampler.IsDone() && uniformSampler.NbPoints() == sampleCount;

            for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
                double parameter = 0.0;
                if (useUniformAbscissa) {
                    const int uniformIndex = startParameter <= endParameter
                        ? sampleIndex + 1
                        : sampleCount - sampleIndex;
                    parameter = uniformSampler.Parameter(uniformIndex);
                }
                else {
                    const double normalized = static_cast<double>(sampleIndex) /
                        static_cast<double>(sampleCount - 1);
                    parameter = startParameter + (endParameter - startParameter) * normalized;
                }
                sampledPoints.push_back(curve.Value(parameter));
            }
            return true;
        };

        for (std::size_t edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
            const double lengthA = GeometryUtils::ComputeEdgeLength(candidate.chainA.edges[edgeIndex]);
            const double lengthB = GeometryUtils::ComputeEdgeLength(candidate.chainB.edges[edgeIndex]);
            if (lengthA <= Precision::Confusion() || lengthB <= Precision::Confusion()) {
                errorMessage = "candidate contains a zero-length source edge pair";
                return false;
            }

            const int sampleCount = computeMidlineSampleCount(lengthA, lengthB);
            std::vector<gp_Pnt> sampledA;
            std::vector<gp_Pnt> sampledB;
            if (!sampleOrientedEdge(
                    candidate.chainA, edgeIndex, sampleCount, sampledA, "chainA") ||
                !sampleOrientedEdge(
                    candidate.chainB, edgeIndex, sampleCount, sampledB, "chainB")) {
                return false;
            }

            std::vector<gp_Pnt> guidancePoints;
            guidancePoints.reserve(static_cast<std::size_t>(sampleCount));
            const gp_Pnt startPoint = BRep_Tool::Pnt(midVertices[edgeIndex]);
            const gp_Pnt endPoint = BRep_Tool::Pnt(midVertices[edgeIndex + 1]);
            guidancePoints.push_back(startPoint);
            for (int sampleIndex = 1; sampleIndex + 1 < sampleCount; ++sampleIndex) {
                const gp_Pnt midpoint(
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].X()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].X()),
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].Y()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].Y()),
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].Z()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].Z()));
                if (guidancePoints.back().Distance(midpoint) > Precision::Confusion() &&
                    midpoint.Distance(endPoint) > Precision::Confusion()) {
                    guidancePoints.push_back(midpoint);
                }
            }
            if (guidancePoints.back().Distance(endPoint) <= Precision::Confusion()) {
                errorMessage = "midline curve guidance points collapse at the segment endpoint";
                return false;
            }
            guidancePoints.push_back(endPoint);

            Handle(TColgp_HArray1OfPnt) interpolationPoints =
                new TColgp_HArray1OfPnt(1, static_cast<Standard_Integer>(guidancePoints.size()));
            for (std::size_t pointIndex = 0; pointIndex < guidancePoints.size(); ++pointIndex) {
                interpolationPoints->SetValue(
                    static_cast<Standard_Integer>(pointIndex + 1),
                    guidancePoints[pointIndex]);
            }

            // 采样中点只约束曲线形状，不创建额外拓扑顶点或 edge。
            GeomAPI_Interpolate interpolator(
                interpolationPoints,
                Standard_False,
                Precision::Confusion());
            interpolator.Perform();
            if (!interpolator.IsDone() || interpolator.Curve().IsNull()) {
                errorMessage = "failed to interpolate a midline curve from sampled midpoints";
                return false;
            }

            Handle(Geom_Curve) midlineCurve = interpolator.Curve();
            BRepBuilderAPI_MakeEdge edgeMaker(
                midlineCurve,
                midVertices[edgeIndex],
                midVertices[edgeIndex + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create midline edge";
                return false;
            }

            const TopoDS_Edge edge = edgeMaker.Edge();
            midlineData.chain.edges.push_back(edge);
            midlineData.chain.length += GeometryUtils::ComputeEdgeLength(edge);
            const std::size_t midlineEdgeIndex = midlineData.chain.edges.size() - 1;

            // source edge 与中线 edge 一一对应，采样密度不再改变拓扑分段数量。
            midlineData.replacementsForChainA[edgeIndex].push_back(midlineEdgeIndex);
            midlineData.replacementsForChainB[edgeIndex].push_back(midlineEdgeIndex);
        }

        if (midlineData.chain.edges.empty() ||
            midlineData.chain.length <= Precision::Confusion()) {
            errorMessage = "midline chain is empty";
            return false;
        }

        for (std::size_t i = 0; i < midlineData.replacementsForChainA.size(); ++i) {
            if (midlineData.replacementsForChainA[i].empty()) {
                errorMessage = "chainA edge has no corresponding midline segment";
                return false;
            }
        }
        for (std::size_t i = 0; i < midlineData.replacementsForChainB.size(); ++i) {
            if (midlineData.replacementsForChainB[i].empty()) {
                errorMessage = "chainB edge has no corresponding midline segment";
                return false;
            }
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "midline chain construction failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::tryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    if (params_.useMidlineCollapse) {
        std::string midlineError;
        if (tryRepairByMidlineCollapse(currentShape, candidate, repairedShape, midlineError)) {
            return true;
        }

        if (!params_.fallbackToSideReplacement) {
            errorMessage = "Closed narrow gap repair failed:\n"
                + RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::StrategyFailureBlock("MidlineCollapse", midlineError), "  ");
            return false;
        }

        std::string sideError;
        if (tryRepairBySideReplacement(currentShape, candidate, repairedShape, sideError)) {
            appendMessage(
                "Warning: repaired by side replacement after earlier strategy failed:\n"
                + RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::StrategyFailureBlock("MidlineCollapse", midlineError), "  "));
            return true;
        }

        errorMessage = "Closed narrow gap repair failed:\n"
            + RepairDiagnostics::IndentMultiline(
                RepairDiagnostics::StrategyFailureBlock("MidlineCollapse", midlineError), "  ")
            + "\n"
            + RepairDiagnostics::IndentMultiline(
                RepairDiagnostics::StrategyFailureBlock("SideReplacement", sideError), "  ");
        return false;
    }

    return tryRepairBySideReplacement(currentShape, candidate, repairedShape, errorMessage);
}

bool HealClosedNarrowGap::tryRepairByMidlineCollapse(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    MidlineChainData midlineData;
    if (!buildMidlineChain(candidate, midlineData, errorMessage)) {
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    if (!rebuildFacesByCollapsingToMidline(
            currentShape,
            candidate,
            midlineData,
            trialShape,
            errorMessage)) {
        return false;
    }

    appendMessage(
        "Midline adjacency after local rebuild: "
        + describeEdgeFaceCounts(trialShape, midlineData.chain.edges));

    if (params_.runShapeFixAfterRepair) {
        try {
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(params_.midlineRepairTolerance);
            fixer.SetMaxTolerance(std::max(params_.midlineRepairTolerance, params_.maxGapWidth));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();

            const TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix_Shape failed after midline collapse: " + failureText(failure);
            return false;
        }
    }

    appendMessage(
        std::string(params_.runShapeFixAfterRepair
            ? "Midline adjacency after final ShapeFix: "
            : "Midline adjacency before validation (final ShapeFix disabled): ")
        + describeEdgeFaceCounts(trialShape, midlineData.chain.edges));

    if (!validateRepairResult(currentShape, trialShape, midlineData.chain, freeEdgesBefore, errorMessage)) {
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by midline collapse: midlineEdges="
           << midlineData.chain.edges.size()
           << ", midlineLength=" << midlineData.chain.length
           << ", repairTolerance=" << params_.midlineRepairTolerance;
        appendMessage(os.str());
    }
    repairedShape = trialShape;
    return true;
}

bool HealClosedNarrowGap::tryRepairBySideReplacement(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    ChainChoice choice = chooseChainToKeep(candidate);
    if (choice.keptChain.edges.empty() || choice.replaceChain.edges.empty()) {
        errorMessage = "invalid chain choice";
        return false;
    }

    if (choice.keptChain.edges.size() != choice.replaceChain.edges.size()) {
        std::ostringstream os;
        os << "first version only supports one-to-one chain edge replacement: keptEdges="
           << choice.keptChain.edges.size()
           << ", replaceEdges=" << choice.replaceChain.edges.size();
        errorMessage = os.str();
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    if (!rebuildFacesByReplacingChain(
            currentShape,
            choice.keptChain,
            choice.replaceChain,
            trialShape,
            errorMessage)) {
        return false;
    }

    if (params_.runShapeFixAfterRepair) {
        try {
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(std::max(Precision::Confusion(), params_.maxGapWidth));
            fixer.SetMaxTolerance(std::max(Precision::Confusion(), params_.maxGapWidth * 10.0));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();

            const TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix_Shape failed after closed gap repair: " + failureText(failure);
            return false;
        }
    }

    if (!validateRepairResult(currentShape, trialShape, choice.keptChain, freeEdgesBefore, errorMessage)) {
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by side replacement: keptEdges="
           << choice.keptChain.edges.size()
           << ", replacedEdges=" << choice.replaceChain.edges.size()
           << ", keptLength=" << choice.keptChain.length
           << ", replacedLength=" << choice.replaceChain.length;
        appendMessage(os.str());
    }
    repairedShape = trialShape;
    return true;
}

HealClosedNarrowGap::ChainChoice HealClosedNarrowGap::chooseChainToKeep(
    const ClosedGapCandidate& candidate) const
{
    ChainChoice choice;
    bool keepA = true;

    if (params_.preferLongerChainToKeep) {
        keepA = candidate.chainA.length >= candidate.chainB.length;
    }
    else if (candidate.chainA.edges.size() != candidate.chainB.edges.size()) {
        keepA = candidate.chainA.edges.size() <= candidate.chainB.edges.size();
    }

    if (keepA) {
        choice.keptChain = candidate.chainA;
        choice.replaceChain = candidate.chainB;
    }
    else {
        choice.keptChain = candidate.chainB;
        choice.replaceChain = candidate.chainA;
    }
    return choice;
}

bool HealClosedNarrowGap::rebuildFacesByCollapsingToMidline(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    const MidlineChainData& midlineData,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }
    if (midlineData.chain.edges.empty()) {
        errorMessage = "midline chain is empty";
        return false;
    }

    std::vector<FaceReplacementPlan> facePlans;
    auto findOrAddFacePlan = [&](const TopoDS_Face& face) -> FaceReplacementPlan* {
        for (FaceReplacementPlan& plan : facePlans) {
            if (plan.face.IsSame(face)) {
                return &plan;
            }
        }

        FaceReplacementPlan plan;
        plan.face = face;
        facePlans.push_back(plan);
        return &facePlans.back();
    };

    // 按原 edge 的 owner face 登记局部替换，不再要求完整 chain 位于单个 face/wire 中。
    auto addChainTasks = [&](const EdgeChain& sourceChain,
                             const std::vector<std::vector<std::size_t>>& replacements,
                             const char* chainName) {
        if (sourceChain.ownerFaces.size() != sourceChain.edges.size() ||
            sourceChain.vertices.size() != sourceChain.edges.size() + 1 ||
            replacements.size() != sourceChain.edges.size()) {
            errorMessage = std::string(chainName) + " replacement mapping is inconsistent";
            return false;
        }

        for (std::size_t edgeIndex = 0; edgeIndex < sourceChain.edges.size(); ++edgeIndex) {
            const TopoDS_Face& ownerFace = sourceChain.ownerFaces[edgeIndex];
            if (ownerFace.IsNull()) {
                errorMessage = std::string(chainName) + " contains an edge without owner face";
                return false;
            }
            if (replacements[edgeIndex].empty()) {
                errorMessage = std::string(chainName) + " edge has no midline replacement segment";
                return false;
            }

            FaceReplacementPlan* plan = findOrAddFacePlan(ownerFace);
            for (const EdgeReplacementTask& existingTask : plan->edgeTasks) {
                if (existingTask.sourceEdge.IsSame(sourceChain.edges[edgeIndex])) {
                    errorMessage = "the same source edge appears more than once in midline replacement tasks";
                    return false;
                }
            }

            EdgeReplacementTask task;
            task.sourceEdge = sourceChain.edges[edgeIndex];
            task.sourceStart = sourceChain.vertices[edgeIndex];
            task.sourceEnd = sourceChain.vertices[edgeIndex + 1];
            task.replacementEdgeIndices = replacements[edgeIndex];
            plan->edgeTasks.push_back(task);
        }
        return true;
    };

    if (!addChainTasks(
            candidate.chainA,
            midlineData.replacementsForChainA,
            "chainA") ||
        !addChainTasks(
            candidate.chainB,
            midlineData.replacementsForChainB,
            "chainB")) {
        return false;
    }

    if (facePlans.empty()) {
        errorMessage = "closed-gap chains have no owner faces";
        return false;
    }

    try {
        BRepTools_ReShape topologyReshaper;
        BRepTools_ReShape vertexReshaper;
        std::vector<TopoDS_Vertex> replacedSourceVertices;
        TopTools_MapOfShape uniqueSourceVertices;

        auto registerVertexReplacement = [&](const TopoDS_Vertex& sourceVertex,
                                             const TopoDS_Vertex& midlineVertex) {
            if (sourceVertex.IsNull() || midlineVertex.IsNull()) {
                errorMessage = "cannot register a null midline vertex replacement";
                return false;
            }
            if (sourceVertex.IsSame(midlineVertex)) {
                return true;
            }

            // ReShape 的替换键和值必须去除 edge use 携带的局部方向，
            // 否则 REVERSED source vertex 会使 replacement use 的方向整体翻转。
            TopoDS_Vertex normalizedSourceVertex = sourceVertex;
            normalizedSourceVertex.Orientation(TopAbs_FORWARD);
            TopoDS_Vertex normalizedMidlineVertex = midlineVertex;
            normalizedMidlineVertex.Orientation(TopAbs_FORWARD);

            if (vertexReshaper.IsRecorded(normalizedSourceVertex)) {
                const TopoDS_Shape existingReplacement =
                    vertexReshaper.Value(normalizedSourceVertex);
                if (existingReplacement.IsNull() ||
                    !existingReplacement.IsSame(normalizedMidlineVertex)) {
                    errorMessage = "one source vertex maps to different midline vertices";
                    return false;
                }
                return true;
            }

            // 两侧对应旧顶点统一替换成同一个规范中点 vertex。
            vertexReshaper.Replace(normalizedSourceVertex, normalizedMidlineVertex);
            topologyReshaper.Replace(normalizedSourceVertex, normalizedMidlineVertex);
            if (uniqueSourceVertices.Add(normalizedSourceVertex)) {
                replacedSourceVertices.push_back(normalizedSourceVertex);
            }
            return true;
        };

        if (candidate.chainA.vertices.size() != midlineData.chain.vertices.size() ||
            candidate.chainB.vertices.size() != midlineData.chain.vertices.size()) {
            errorMessage = "midline vertex mapping is inconsistent with candidate chains";
            return false;
        }
        for (std::size_t vertexIndex = 1;
             vertexIndex + 1 < midlineData.chain.vertices.size();
             ++vertexIndex) {
            if (!registerVertexReplacement(
                    candidate.chainA.vertices[vertexIndex],
                    midlineData.chain.vertices[vertexIndex]) ||
                !registerVertexReplacement(
                    candidate.chainB.vertices[vertexIndex],
                    midlineData.chain.vertices[vertexIndex])) {
                return false;
            }
        }

        // 先为所有连接旧内部顶点的 edge 创建唯一替代品，确保相邻 face 继续共享同一 edge TShape。
        TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap;
        TopologyUtils::BuildVertexEdgeMap(currentShape, vertexEdgeMap);
        TopTools_MapOfShape transformedEdges;
        for (const TopoDS_Vertex& sourceVertex : replacedSourceVertices) {
            if (!vertexEdgeMap.Contains(sourceVertex)) {
                continue;
            }

            const TopTools_ListOfShape& incidentEdges = vertexEdgeMap.FindFromKey(sourceVertex);
            for (TopTools_ListIteratorOfListOfShape edgeIt(incidentEdges);
                 edgeIt.More();
                 edgeIt.Next()) {
                const TopoDS_Edge sourceEdge = TopoDS::Edge(edgeIt.Value());
                if (sourceEdge.IsNull() || !transformedEdges.Add(sourceEdge)) {
                    continue;
                }

                const TopoDS_Shape transformedShape = vertexReshaper.Apply(sourceEdge);
                if (transformedShape.IsNull() || transformedShape.ShapeType() != TopAbs_EDGE) {
                    errorMessage = "failed to propagate a midpoint vertex to an incident edge";
                    return false;
                }
                TopoDS_Edge transformedEdge = TopoDS::Edge(transformedShape);
                transformedEdge.Orientation(sourceEdge.Orientation());
                topologyReshaper.Replace(sourceEdge, transformedEdge);
            }
        }

        const double tolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
        for (const FaceReplacementPlan& plan : facePlans) {
            TopoDS_Face rebuiltFace;
            if (!rebuildFaceByReplacingMappedEdges(
                    plan.face,
                    plan.edgeTasks,
                    midlineData.chain,
                    topologyReshaper,
                    tolerance,
                    rebuiltFace,
                    errorMessage)) {
                return false;
            }

            topologyReshaper.Replace(plan.face, rebuiltFace);
        }

        repairedShape = topologyReshaper.Apply(currentShape);
        if (repairedShape.IsNull()) {
            errorMessage = "BRepTools_ReShape produced a null shape after midline collapse";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "midline shape replacement failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFaceByReplacingMappedEdges(
    const TopoDS_Face& faceToRebuild,
    const std::vector<EdgeReplacementTask>& edgeTasks,
    const EdgeChain& midlineChain,
    BRepTools_ReShape& topologyReshaper,
    double tolerance,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage)
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull()) {
        errorMessage = "face to rebuild is null";
        return false;
    }
    if (edgeTasks.empty()) {
        errorMessage = "face has no midline edge replacement tasks";
        return false;
    }
    if (midlineChain.edges.empty()) {
        errorMessage = "canonical midline chain is empty";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "face to rebuild has no support surface";
            return false;
        }

        const double faceTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(faceToRebuild),
            tolerance
        });

        BRep_Builder builder;
        TopoDS_Face rebuilt;
        builder.MakeFace(rebuilt, surface, surfaceLocation, faceTolerance);
        builder.NaturalRestriction(rebuilt, BRep_Tool::NaturalRestriction(forwardFace));
        rebuilt.Orientation(TopAbs_FORWARD);

        std::vector<bool> replacedTasks(edgeTasks.size(), false);
        auto findTaskIndex = [&](const TopoDS_Edge& edge) {
            for (std::size_t taskIndex = 0; taskIndex < edgeTasks.size(); ++taskIndex) {
                if (edgeTasks[taskIndex].sourceEdge.IsSame(edge)) {
                    return static_cast<int>(taskIndex);
                }
            }
            return -1;
        };

        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            // 直接子 edge use 是旧 wire 的完整基准，不能只依赖 WireExplorer 的返回结果。
            std::vector<TopoDS_Edge> directEdgeUses;
            bool directWireHasReplacement = false;
            for (TopoDS_Iterator directIterator(oldWire, Standard_False, Standard_False);
                 directIterator.More();
                 directIterator.Next()) {
                if (directIterator.Value().ShapeType() != TopAbs_EDGE) {
                    continue;
                }

                const TopoDS_Edge directEdge = TopoDS::Edge(directIterator.Value());
                directEdgeUses.push_back(directEdge);
                directWireHasReplacement =
                    directWireHasReplacement || findTaskIndex(directEdge) >= 0;
            }

            std::vector<TopoDS_Edge> wireEdges;
            bool wireHasReplacement = false;
            for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                const TopoDS_Edge edgeInWire = edgeExplorer.Current();
                wireEdges.push_back(edgeInWire);
                wireHasReplacement = wireHasReplacement || findTaskIndex(edgeInWire) >= 0;
            }

            if (directWireHasReplacement) {
                // 按 IsSame 做多重集匹配，既检查数量，也检查 WireExplorer 是否遗漏具体 edge use。
                std::vector<bool> matchedDirectUses(directEdgeUses.size(), false);
                int unexpectedExploredUses = 0;
                for (const TopoDS_Edge& exploredEdge : wireEdges) {
                    bool matched = false;
                    for (std::size_t directIndex = 0;
                         directIndex < directEdgeUses.size();
                         ++directIndex) {
                        if (!matchedDirectUses[directIndex] &&
                            directEdgeUses[directIndex].IsSame(exploredEdge)) {
                            matchedDirectUses[directIndex] = true;
                            matched = true;
                            break;
                        }
                    }
                    if (!matched) {
                        ++unexpectedExploredUses;
                    }
                }

                const int missingDirectUses = static_cast<int>(std::count(
                    matchedDirectUses.begin(), matchedDirectUses.end(), false));
                {
                    std::ostringstream os;
                    os << "Midline wire traversal check: directEdgeUses="
                       << directEdgeUses.size()
                       << ", exploredEdgeUses=" << wireEdges.size()
                       << ", missingDirectUses=" << missingDirectUses
                       << ", unexpectedExploredUses=" << unexpectedExploredUses;
                    appendMessage(os.str());
                }

                if (missingDirectUses > 0 || unexpectedExploredUses > 0 ||
                    !wireHasReplacement) {
                    std::ostringstream os;
                    os << "BRepTools_WireExplorer did not cover every direct edge use of the affected wire: "
                       << "direct=" << directEdgeUses.size()
                       << ", explored=" << wireEdges.size()
                       << ", missing=" << missingDirectUses
                       << ", unexpected=" << unexpectedExploredUses;
                    errorMessage = os.str();
                    return false;
                }
            }

            TopoDS_Wire wireToAdd = oldWire;
            if (wireHasReplacement) {
                BRepBuilderAPI_MakeWire wireMaker;
                for (const TopoDS_Edge& edgeInWire : wireEdges) {
                    const int taskIndex = findTaskIndex(edgeInWire);
                    if (taskIndex < 0) {
                        const TopoDS_Shape adjustedShape = topologyReshaper.Apply(edgeInWire);
                        if (adjustedShape.IsNull() || adjustedShape.ShapeType() != TopAbs_EDGE) {
                            errorMessage = "failed to apply midpoint vertices to a neighboring wire edge";
                            return false;
                        }
                        TopoDS_Edge adjustedEdge = TopoDS::Edge(adjustedShape);
                        adjustedEdge.Orientation(edgeInWire.Orientation());
                        wireMaker.Add(adjustedEdge);
                        continue;
                    }

                    const std::size_t mappedIndex = static_cast<std::size_t>(taskIndex);
                    if (replacedTasks[mappedIndex]) {
                        errorMessage = "a source edge appears more than once in face wires";
                        return false;
                    }

                    const EdgeReplacementTask& task = edgeTasks[mappedIndex];
                    TopoDS_Vertex wireStart;
                    TopoDS_Vertex wireEnd;
                    if (!getEdgeVertices(edgeInWire, wireStart, wireEnd)) {
                        errorMessage = "failed to obtain source edge vertices from face wire";
                        return false;
                    }

                    bool reverseReplacement = false;
                    if (wireStart.IsSame(task.sourceStart) && wireEnd.IsSame(task.sourceEnd)) {
                        reverseReplacement = false;
                    }
                    else if (wireStart.IsSame(task.sourceEnd) && wireEnd.IsSame(task.sourceStart)) {
                        reverseReplacement = true;
                    }
                    else {
                        errorMessage = "source edge vertices do not match its chain interval";
                        return false;
                    }

                    auto addReplacementUse = [&](std::size_t replacementEdgeIndex,
                                                 bool reverseUse) {
                        if (replacementEdgeIndex >= midlineChain.edges.size()) {
                            errorMessage = "mapped midline edge index is out of range";
                            return false;
                        }

                        // 所有 face 都从规范化中线 chain 取同一个底层 edge，只调整当前 wire use 的方向。
                        TopoDS_Edge replacementUse =
                            midlineChain.edges[replacementEdgeIndex];
                        if (reverseUse) {
                            replacementUse.Reverse();
                        }

                        // 先为当前 face 补建并验证 pcurve，再把同一个底层中线 edge use 加入 wire。
                        if (!ensurePCurveOnFace(
                                replacementUse,
                                forwardFace,
                                faceTolerance,
                                errorMessage)) {
                            errorMessage = "failed to repair mapped midline pcurve: " + errorMessage;
                            return false;
                        }
                        wireMaker.Add(replacementUse);
                        return true;
                    };

                    if (reverseReplacement) {
                        for (int replacementIndex =
                                 static_cast<int>(task.replacementEdgeIndices.size()) - 1;
                             replacementIndex >= 0;
                             --replacementIndex) {
                            if (!addReplacementUse(
                                    task.replacementEdgeIndices[
                                        static_cast<std::size_t>(replacementIndex)],
                                    true)) {
                                return false;
                            }
                        }
                    }
                    else {
                        for (const std::size_t replacementEdgeIndex :
                             task.replacementEdgeIndices) {
                            if (!addReplacementUse(replacementEdgeIndex, false)) {
                                return false;
                            }
                        }
                    }

                    replacedTasks[mappedIndex] = true;
                }

                if (!wireMaker.IsDone()) {
                    errorMessage = "failed to rebuild wire while replacing mapped midline edges";
                    return false;
                }

                wireToAdd = wireMaker.Wire();
                wireToAdd.Orientation(oldWire.Orientation());

                // 局部修复仅处理受影响 wire，最终仍由 ShapeFix_Shape 统一检查和修复整体 shape。
                ShapeFix_Wire wireFix(wireToAdd, forwardFace, faceTolerance);
                wireFix.ModifyTopologyMode() = Standard_False;
                wireFix.ModifyGeometryMode() = Standard_True;
                wireFix.ClosedWireMode() = Standard_True;
                wireFix.FixReorderMode() = 1;
                wireFix.FixConnectedMode() = 1;
                wireFix.FixEdgeCurvesMode() = 1;
                wireFix.FixReversed2dMode() = 1;
                wireFix.FixSameParameterMode() = 1;
                wireFix.Perform();

                TopoDS_Wire fixedWire = wireFix.WireAPIMake();
                if (fixedWire.IsNull()) {
                    fixedWire = wireFix.Wire();
                }
                if (fixedWire.IsNull()) {
                    errorMessage = "ShapeFix_Wire returned null wire after mapped midline replacement";
                    return false;
                }

                fixedWire.Orientation(oldWire.Orientation());
                wireToAdd = fixedWire;
            }
            else {
                const TopoDS_Shape adjustedWireShape = topologyReshaper.Apply(oldWire);
                if (adjustedWireShape.IsNull() || adjustedWireShape.ShapeType() != TopAbs_WIRE) {
                    errorMessage = "failed to propagate midpoint vertices to an unchanged face wire";
                    return false;
                }
                wireToAdd = TopoDS::Wire(adjustedWireShape);
                wireToAdd.Orientation(oldWire.Orientation());
            }

            builder.Add(rebuilt, wireToAdd);
        }

        for (std::size_t taskIndex = 0; taskIndex < replacedTasks.size(); ++taskIndex) {
            if (!replacedTasks[taskIndex]) {
                errorMessage = "face does not contain one of its mapped source edges";
                return false;
            }
        }

        // 恢复 face 级局部修复，修正 wire 方向以及新中线 edge use 的曲面边界表示。
        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTolerance);
        faceFix.SetMaxTolerance(std::max(faceTolerance, params_.maxGapWidth));
        faceFix.FixWireMode() = 1;
        faceFix.FixOrientationMode() = 1;
        faceFix.FixMissingSeamMode() = 1;
        faceFix.FixAddNaturalBoundMode() = 0;
        faceFix.Perform();

        const TopoDS_Shape fixedShape = faceFix.Result();
        if (!fixedShape.IsNull() && fixedShape.ShapeType() == TopAbs_FACE) {
            rebuiltFace = TopoDS::Face(fixedShape);
        }
        else {
            rebuiltFace = faceFix.Face();
        }
        if (rebuiltFace.IsNull()) {
            errorMessage = "ShapeFix_Face returned null face after mapped midline replacement";
            return false;
        }

        rebuiltFace.Orientation(faceToRebuild.Orientation());
    }
    catch (Standard_Failure& failure) {
        errorMessage = "face rebuild with mapped midline failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFacesByReplacingChain(
    const TopoDS_Shape& currentShape,
    const EdgeChain& keptChain,
    const EdgeChain& replaceChain,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    std::vector<TopoDS_Face> facesToRebuild;
    for (const TopoDS_Face& ownerFace : replaceChain.ownerFaces) {
        if (!ownerFace.IsNull() && !isSameFaceInVector(facesToRebuild, ownerFace)) {
            facesToRebuild.push_back(ownerFace);
        }
    }

    if (facesToRebuild.empty()) {
        errorMessage = "replace chain has no owner faces";
        return false;
    }

    try {
        BRepTools_ReShape reshaper;
        const double tolerance = std::max(params_.maxGapWidth, Precision::Confusion());
        for (const TopoDS_Face& faceToRebuild : facesToRebuild) {
            TopoDS_Face rebuiltFace;
            if (!rebuildFaceByReplacingChainEdges(
                    faceToRebuild,
                    keptChain,
                    replaceChain,
                    tolerance,
                    rebuiltFace,
                    errorMessage)) {
                return false;
            }

            reshaper.Replace(faceToRebuild, rebuiltFace);
        }

        repairedShape = reshaper.Apply(currentShape);
        if (repairedShape.IsNull()) {
            errorMessage = "BRepTools_ReShape produced a null shape";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "shape replacement failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFaceByReplacingChainEdges(
    const TopoDS_Face& faceToRebuild,
    const EdgeChain& keptChain,
    const EdgeChain& replaceChain,
    double tolerance,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage) const
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull()) {
        errorMessage = "face to rebuild is null";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "face to rebuild has no support surface";
            return false;
        }

        const double faceTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(faceToRebuild),
            tolerance
        });

        BRep_Builder builder;
        TopoDS_Face rebuilt;
        builder.MakeFace(rebuilt, surface, surfaceLocation, faceTolerance);
        builder.NaturalRestriction(rebuilt, BRep_Tool::NaturalRestriction(forwardFace));
        rebuilt.Orientation(TopAbs_FORWARD);

        bool replacedAnyEdge = false;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            bool wireContainsReplaceEdge = false;
            for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                if (findSameEdgeIndex(replaceChain.edges, edgeExplorer.Current()) >= 0) {
                    wireContainsReplaceEdge = true;
                    break;
                }
            }

            TopoDS_Wire wireToAdd = oldWire;
            if (wireContainsReplaceEdge) {
                BRepBuilderAPI_MakeWire wireMaker;
                for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                     edgeExplorer.More();
                     edgeExplorer.Next()) {
                    const TopoDS_Edge edgeInWire = edgeExplorer.Current();
                    const int replaceIndex = findSameEdgeIndex(replaceChain.edges, edgeInWire);
                    if (replaceIndex < 0) {
                        wireMaker.Add(edgeInWire);
                        continue;
                    }

                    TopoDS_Edge replacementUse = keptChain.edges[static_cast<std::size_t>(replaceIndex)];
                    replacementUse.Orientation(edgeInWire.Orientation());
                    if (!ensurePCurveOnFace(replacementUse, forwardFace, faceTolerance, errorMessage)) {
                        errorMessage = "failed to repair replacement pcurve: " + errorMessage;
                        return false;
                    }

                    wireMaker.Add(replacementUse);
                    replacedAnyEdge = true;
                }

                if (!wireMaker.IsDone()) {
                    errorMessage = "failed to rebuild wire while replacing chain edges";
                    return false;
                }

                wireToAdd = wireMaker.Wire();
                wireToAdd.Orientation(oldWire.Orientation());

                ShapeFix_Wire wireFix(wireToAdd, forwardFace, faceTolerance);
                wireFix.ModifyTopologyMode() = Standard_False;
                wireFix.ModifyGeometryMode() = Standard_True;
                wireFix.ClosedWireMode() = Standard_True;
                wireFix.FixReorderMode() = 1;
                wireFix.FixConnectedMode() = 1;
                wireFix.FixEdgeCurvesMode() = 1;
                wireFix.FixReversed2dMode() = 1;
                wireFix.FixSameParameterMode() = 1;
                wireFix.Perform();

                TopoDS_Wire fixedWire = wireFix.WireAPIMake();
                if (fixedWire.IsNull()) {
                    fixedWire = wireFix.Wire();
                }
                if (fixedWire.IsNull()) {
                    errorMessage = "ShapeFix_Wire returned null wire";
                    return false;
                }

                fixedWire.Orientation(oldWire.Orientation());
                wireToAdd = fixedWire;
            }

            builder.Add(rebuilt, wireToAdd);
        }

        if (!replacedAnyEdge) {
            errorMessage = "face does not contain replace chain edges";
            return false;
        }

        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTolerance);
        faceFix.SetMaxTolerance(std::max(faceTolerance, params_.maxGapWidth * 10.0));
        faceFix.FixWireMode() = 1;
        faceFix.FixOrientationMode() = 1;
        faceFix.FixMissingSeamMode() = 1;
        faceFix.FixAddNaturalBoundMode() = 0;
        faceFix.Perform();

        const TopoDS_Shape fixedShape = faceFix.Result();
        if (!fixedShape.IsNull() && fixedShape.ShapeType() == TopAbs_FACE) {
            rebuiltFace = TopoDS::Face(fixedShape);
        }
        else {
            rebuiltFace = faceFix.Face();
        }

        if (rebuiltFace.IsNull()) {
            errorMessage = "ShapeFix_Face returned null face";
            return false;
        }

        rebuiltFace.Orientation(faceToRebuild.Orientation());
    }
    catch (Standard_Failure& failure) {
        errorMessage = "face rebuild failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::validateRepairResult(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& afterShape,
    const EdgeChain& keptChain,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (beforeShape.IsNull() || afterShape.IsNull()) {
        errorMessage = "before or after shape is null";
        return false;
    }

    std::string existingError;
    const RepairDiagnostics::FailureList existingFailures =
        RepairDiagnostics::CollectBRepCheckFailures(beforeShape, &existingError);
    if (!existingError.empty()) {
        errorMessage = "existing BRepCheck failure collection failed: " + existingError;
        return false;
    }

    std::string repairedError;
    const RepairDiagnostics::FailureList repairedFailures =
        RepairDiagnostics::CollectBRepCheckFailures(afterShape, &repairedError);
    if (!repairedError.empty()) {
        errorMessage = "repaired BRepCheck failure collection failed: " + repairedError;
        return false;
    }

    // face 重建会改变 TShape 身份；按错误类型和状态的数量比较，避免把旧错误迁移误判为新增。
    const RepairDiagnostics::FailureList newFailures =
        RepairDiagnostics::FindBRepCheckFailureCountIncrease(
            existingFailures,
            repairedFailures);

    if (!newFailures.empty() && params_.rejectInvalidResult) {
        errorMessage = "new BRepCheck failures after closed narrow gap repair:"
            + RepairDiagnostics::DescribeBRepCheckFailures(newFailures);
        return false;
    }

    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (freeEdgesAfter > freeEdgesBefore && params_.rejectFreeEdgeIncrease) {
        std::ostringstream os;
        os << "free edge count increased after closed narrow gap repair: before="
           << freeEdgesBefore << ", after=" << freeEdgesAfter;
        errorMessage = os.str();
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(afterShape, edgeFaceMap);
    for (const TopoDS_Edge& keptEdge : keptChain.edges) {
        int faceCount = 0;
        for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
            const TopoDS_Edge mappedEdge = TopoDS::Edge(edgeFaceMap.FindKey(i));
            if (!mappedEdge.IsNull() && mappedEdge.IsSame(keptEdge)) {
                faceCount += listExtent(edgeFaceMap.FindFromIndex(i));
            }
        }

        if (faceCount < 2) {
            errorMessage = "kept chain edge is not shared after repair";
            return false;
        }
    }

    return true;
}

int HealClosedNarrowGap::computeAutomaticIterationLimit(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return 1;
    }

    TopTools_IndexedMapOfShape edgeMap;
    TopExp::MapShapes(shape, TopAbs_EDGE, edgeMap);
    if (edgeMap.Extent() <= 0) {
        return 1;
    }

    return std::min(edgeMap.Extent(), kMaxAutomaticRepairIterations);
}

int HealClosedNarrowGap::computeMidlineSampleCount(double lengthA, double lengthB) const
{
    const int fallbackCount = std::max(2, params_.sampleCount);
    if (!params_.useAdaptiveSampleCount ||
        params_.targetSampleSpacing <= Precision::Confusion()) {
        return fallbackCount;
    }

    const double pairLength = std::max(lengthA, lengthB);
    if (!std::isfinite(pairLength) || pairLength <= Precision::Confusion()) {
        return fallbackCount;
    }

    // 加一是为了让采样序列同时包含 source edge 的两个端点。
    const int countByLength = static_cast<int>(std::ceil(
        pairLength / params_.targetSampleSpacing)) + 1;
    return std::max(
        params_.minSampleCount,
        std::min(params_.maxSampleCount, countByLength));
}

int HealClosedNarrowGap::findSameEdgeIndex(
    const std::vector<TopoDS_Edge>& edges,
    const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return -1;
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!edges[i].IsNull() && edges[i].IsSame(edge)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool HealClosedNarrowGap::ensurePCurveOnFace(
    TopoDS_Edge& edge,
    const TopoDS_Face& face,
    double tolerance,
    std::string& errorMessage) const
{
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "edge or face is null";
        return false;
    }

    if (!GeometryUtils::RepairPCurveForEdgeOnFace(
            edge,
            face,
            std::max(tolerance, Precision::Confusion()),
            &errorMessage)) {
        return false;
    }
    return true;
}

void HealClosedNarrowGap::appendMessage(const std::string& message)
{
    messages_.push_back(message);
}
