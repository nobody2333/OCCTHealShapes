#include "HealClosedNarrowGap.h"

#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
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
#include <TopoDS_Wire.hxx>
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

struct BRepCheckFailureInfo {
    TopoDS_Shape shape;
    TopAbs_ShapeEnum type = TopAbs_SHAPE;
    std::vector<BRepCheck_Status> statuses;
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

bool containsStatus(const std::vector<BRepCheck_Status>& statuses, BRepCheck_Status status)
{
    return std::find(statuses.begin(), statuses.end(), status) != statuses.end();
}

bool statusesAreSubset(
    const std::vector<BRepCheck_Status>& statuses,
    const std::vector<BRepCheck_Status>& existingStatuses)
{
    for (const BRepCheck_Status status : statuses) {
        if (!containsStatus(existingStatuses, status)) {
            return false;
        }
    }
    return true;
}

std::vector<BRepCheckFailureInfo> collectBRepCheckFailures(
    const TopoDS_Shape& shape,
    std::string* errorMessage = nullptr)
{
    std::vector<BRepCheckFailureInfo> failures;
    if (shape.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "shape is null";
        }
        return failures;
    }

    try {
        BRepCheck_Analyzer analyzer(shape, Standard_True);
        TopTools_IndexedMapOfShape allShapes;
        TopExp::MapShapes(shape, allShapes);

        for (int i = 1; i <= allShapes.Extent(); ++i) {
            const TopoDS_Shape& subShape = allShapes.FindKey(i);
            Handle(BRepCheck_Result) result = analyzer.Result(subShape);
            if (result.IsNull()) {
                continue;
            }

            BRepCheckFailureInfo failure;
            failure.shape = subShape;
            failure.type = subShape.ShapeType();

            const BRepCheck_ListOfStatus& statuses = result->Status();
            for (BRepCheck_ListIteratorOfListOfStatus it(statuses); it.More(); it.Next()) {
                const BRepCheck_Status status = it.Value();
                if (status == BRepCheck_NoError || containsStatus(failure.statuses, status)) {
                    continue;
                }
                failure.statuses.push_back(status);
            }

            if (!failure.statuses.empty()) {
                failures.push_back(failure);
            }
        }
    }
    catch (Standard_Failure& failure) {
        if (errorMessage != nullptr) {
            *errorMessage = "BRepCheck failure collection failed: " + failureText(failure);
        }
    }

    return failures;
}

bool isExistingBRepCheckFailure(
    const BRepCheckFailureInfo& failure,
    const std::vector<BRepCheckFailureInfo>& existingFailures)
{
    for (const BRepCheckFailureInfo& existingFailure : existingFailures) {
        if (failure.type != existingFailure.type || !failure.shape.IsSame(existingFailure.shape)) {
            continue;
        }

        if (statusesAreSubset(failure.statuses, existingFailure.statuses)) {
            return true;
        }
    }
    return false;
}

std::string describeBRepCheckFailures(const std::vector<BRepCheckFailureInfo>& failures)
{
    std::ostringstream os;
    int reported = 0;
    constexpr int maxReported = 30;
    for (const BRepCheckFailureInfo& failure : failures) {
        if (reported >= maxReported) {
            os << " ...";
            break;
        }

        os << " [type=" << TopologyUtils::ShapeTypeName(failure.type) << " status=";
        for (std::size_t i = 0; i < failure.statuses.size(); ++i) {
            if (i > 0) {
                os << ",";
            }
            os << static_cast<int>(failure.statuses[i]);
        }
        os << "]";
        ++reported;
    }
    return os.str();
}

int listExtent(const TopTools_ListOfShape& shapes)
{
    int count = 0;
    for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
        ++count;
    }
    return count;
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
        os << "Start HealClosedNarrowGap: iterationLimit=" << iterationLimit
           << ", maxGapWidth=" << params_.maxGapWidth
           << ", minLengthRatio=" << params_.minLengthRatio
           << ", minLengthWidthRatio=" << params_.minLengthWidthRatio
           << ", maxTipAngle=" << params_.maxTipAngle
           << ", useMidlineCollapse=" << (params_.useMidlineCollapse ? "true" : "false")
           << ", fallbackToSideReplacement=" << (params_.fallbackToSideReplacement ? "true" : "false")
           << ", midlineRepairTolerance=" << params_.midlineRepairTolerance;
        appendMessage(os.str());
    }

    for (int iter = 0; iter < iterationLimit; ++iter) {
        std::vector<ClosedGapCandidate> candidates = extractCandidates(currentShape);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const ClosedGapCandidate& lhs, const ClosedGapCandidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "Iteration " << (iter + 1)
               << ": closed narrow gap candidate count=" << candidates.size();
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
                appendMessage("Candidate failed: " + errorMessage);
                continue;
            }

            currentShape = repairedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++repairedCount_;
            repairedThisIteration = true;

            std::ostringstream os;
            os << "Repaired closed narrow gap: chainAEdges=" << candidate.chainA.edges.size()
               << ", chainBEdges=" << candidate.chainB.edges.size()
               << ", maxWidth=" << candidate.maxWidth
               << ", meanWidth=" << candidate.meanWidth
               << ", confidence=" << candidate.confidence;
            appendMessage(os.str());
            break;
        }

        if (!repairedThisIteration) {
            break;
        }
    }

    resultShape_ = currentShape;
    {
        std::ostringstream os;
        os << "Finish HealClosedNarrowGap: repairedCount=" << repairedCount_
           << ", modified=" << (modified_ ? "true" : "false");
        appendMessage(os.str());
    }
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
    EdgeChain& midlineChain,
    std::string& errorMessage) const
{
    midlineChain = EdgeChain();
    if (candidate.chainA.edges.empty() || candidate.chainB.edges.empty()) {
        errorMessage = "candidate chains are empty";
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
        // 中线第一版采用折线近似：按归一化长度采样两侧 chain，再取空间中点。
        const int sampleCount = std::max(2, params_.sampleCount);
        std::vector<TopoDS_Vertex> midVertices;
        midVertices.reserve(static_cast<std::size_t>(sampleCount));

        BRep_Builder builder;
        for (int sample = 0; sample < sampleCount; ++sample) {
            const double t = sampleCount == 1
                ? 0.0
                : static_cast<double>(sample) / static_cast<double>(sampleCount - 1);

            gp_Pnt pointA;
            gp_Pnt pointB;
            if (!sampleChainPoint(candidate.chainA, t, pointA) ||
                !sampleChainPoint(candidate.chainB, t, pointB)) {
                errorMessage = "failed to sample candidate chains for midline";
                return false;
            }

            if (sample == 0) {
                // 起点沿用原 tip vertex，保证中线 chain 能和原 wire 的相邻边拓扑相接。
                midVertices.push_back(candidate.chainA.vertices.front());
            }
            else if (sample + 1 == sampleCount) {
                // 终点同样沿用原 tip vertex，避免只靠几何重合连接 wire。
                midVertices.push_back(candidate.chainA.vertices.back());
            }
            else {
                const gp_Pnt midpoint(
                    0.5 * (pointA.X() + pointB.X()),
                    0.5 * (pointA.Y() + pointB.Y()),
                    0.5 * (pointA.Z() + pointB.Z()));

                TopoDS_Vertex midVertex;
                builder.MakeVertex(midVertex, midpoint, params_.midlineRepairTolerance);
                if (midVertex.IsNull()) {
                    errorMessage = "failed to create midline vertex";
                    return false;
                }
                midVertices.push_back(midVertex);
            }
        }

        midlineChain.vertices = midVertices;
        for (std::size_t i = 0; i + 1 < midVertices.size(); ++i) {
            const gp_Pnt p0 = BRep_Tool::Pnt(midVertices[i]);
            const gp_Pnt p1 = BRep_Tool::Pnt(midVertices[i + 1]);
            if (p0.Distance(p1) <= Precision::Confusion()) {
                errorMessage = "midline contains a zero-length segment";
                return false;
            }

            BRepBuilderAPI_MakeEdge edgeMaker(midVertices[i], midVertices[i + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create midline edge";
                return false;
            }

            const TopoDS_Edge edge = edgeMaker.Edge();
            midlineChain.edges.push_back(edge);
            midlineChain.length += GeometryUtils::ComputeEdgeLength(edge);
        }

        if (midlineChain.edges.empty() || midlineChain.length <= Precision::Confusion()) {
            errorMessage = "midline chain is empty";
            return false;
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
            errorMessage = "midline collapse failed: " + midlineError;
            return false;
        }

        appendMessage("Warning: midline collapse failed, fallback to side replacement: " + midlineError);
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

    EdgeChain midlineChain;
    if (!buildMidlineChain(candidate, midlineChain, errorMessage)) {
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    if (!rebuildFacesByCollapsingToMidline(
            currentShape,
            candidate,
            midlineChain,
            trialShape,
            errorMessage)) {
        return false;
    }

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

    if (!validateRepairResult(currentShape, trialShape, midlineChain, freeEdgesBefore, errorMessage)) {
        return false;
    }

    appendMessage("Candidate repaired by midline collapse");
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
    const EdgeChain& midlineChain,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }
    if (midlineChain.edges.empty()) {
        errorMessage = "midline chain is empty";
        return false;
    }

    std::vector<TopoDS_Face> facesToRebuild;
    std::vector<const EdgeChain*> sourceChains;
    auto addFace = [&](const TopoDS_Face& face, const EdgeChain& sourceChain) {
        if (face.IsNull()) {
            return true;
        }

        for (std::size_t i = 0; i < facesToRebuild.size(); ++i) {
            if (!facesToRebuild[i].IsSame(face)) {
                continue;
            }

            if (sourceChains[i] != &sourceChain) {
                errorMessage = "one face contains both closed-gap chains; midline collapse does not support this case";
                return false;
            }
            return true;
        }

        facesToRebuild.push_back(face);
        sourceChains.push_back(&sourceChain);
        return true;
    };

    for (const TopoDS_Face& face : candidate.chainA.ownerFaces) {
        if (!addFace(face, candidate.chainA)) {
            return false;
        }
    }
    for (const TopoDS_Face& face : candidate.chainB.ownerFaces) {
        if (!addFace(face, candidate.chainB)) {
            return false;
        }
    }

    if (facesToRebuild.empty()) {
        errorMessage = "closed-gap chains have no owner faces";
        return false;
    }

    try {
        BRepTools_ReShape reshaper;
        const double tolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
        for (std::size_t i = 0; i < facesToRebuild.size(); ++i) {
            TopoDS_Face rebuiltFace;
            if (!rebuildFaceByReplacingChainWithChain(
                    facesToRebuild[i],
                    *sourceChains[i],
                    midlineChain,
                    tolerance,
                    rebuiltFace,
                    errorMessage)) {
                return false;
            }

            reshaper.Replace(facesToRebuild[i], rebuiltFace);
        }

        repairedShape = reshaper.Apply(currentShape);
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

bool HealClosedNarrowGap::rebuildFaceByReplacingChainWithChain(
    const TopoDS_Face& faceToRebuild,
    const EdgeChain& sourceChain,
    const EdgeChain& replacementChain,
    double tolerance,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage) const
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull()) {
        errorMessage = "face to rebuild is null";
        return false;
    }
    if (sourceChain.edges.empty() || replacementChain.edges.empty()) {
        errorMessage = "source or replacement chain is empty";
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

        bool replacedChain = false;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            std::vector<TopoDS_Edge> wireEdges;
            for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                wireEdges.push_back(edgeExplorer.Current());
            }

            int sourceStartIndex = -1;
            bool sourceReversed = false;
            const bool wireContainsSourceChain =
                findContiguousChainInWire(wireEdges, sourceChain, sourceStartIndex, sourceReversed);

            TopoDS_Wire wireToAdd = oldWire;
            if (wireContainsSourceChain) {
                BRepBuilderAPI_MakeWire wireMaker;
                const int wireEdgeCount = static_cast<int>(wireEdges.size());
                const int sourceEdgeCount = static_cast<int>(sourceChain.edges.size());
                const bool sourceWrapsWireStart = sourceStartIndex + sourceEdgeCount > wireEdgeCount;

                // 按 source chain 在 wire 中的方向加入中线 replacement chain。
                const auto addReplacementChain = [&]() {
                    if (sourceReversed) {
                        for (int replaceIndex = static_cast<int>(replacementChain.edges.size()) - 1;
                             replaceIndex >= 0;
                             --replaceIndex) {
                            TopoDS_Edge replacementUse =
                                replacementChain.edges[static_cast<std::size_t>(replaceIndex)];
                            replacementUse.Reverse();
                            if (!ensurePCurveOnFace(replacementUse, forwardFace, faceTolerance, errorMessage)) {
                                errorMessage = "failed to repair reversed midline pcurve: " + errorMessage;
                                return false;
                            }
                            wireMaker.Add(replacementUse);
                        }
                    }
                    else {
                        for (const TopoDS_Edge& replacementEdge : replacementChain.edges) {
                            TopoDS_Edge replacementUse = replacementEdge;
                            if (!ensurePCurveOnFace(replacementUse, forwardFace, faceTolerance, errorMessage)) {
                                errorMessage = "failed to repair midline pcurve: " + errorMessage;
                                return false;
                            }
                            wireMaker.Add(replacementUse);
                        }
                    }

                    return true;
                };

                if (sourceWrapsWireStart) {
                    // source chain 跨越 BRepTools_WireExplorer 的起点时，闭合 wire 可以从 sourceStartIndex
                    // 重新开始构造：先加入中线 chain，再按环形顺序加入剩余旧 edge。
                    if (!addReplacementChain()) {
                        return false;
                    }

                    for (int offset = sourceEdgeCount; offset < wireEdgeCount; ++offset) {
                        const int edgeIndex = (sourceStartIndex + offset) % wireEdgeCount;
                        wireMaker.Add(wireEdges[static_cast<std::size_t>(edgeIndex)]);
                    }

                    replacedChain = true;
                }
                else {
                    for (int edgeIndex = 0; edgeIndex < wireEdgeCount;) {
                        if (edgeIndex == sourceStartIndex) {
                            // 用同一组中线 edge 替换原 wire 中连续的 source chain。
                            if (!addReplacementChain()) {
                                return false;
                            }

                            edgeIndex += sourceEdgeCount;
                            replacedChain = true;
                            continue;
                        }

                        wireMaker.Add(wireEdges[static_cast<std::size_t>(edgeIndex)]);
                        ++edgeIndex;
                    }
                }

                if (!wireMaker.IsDone()) {
                    errorMessage = "failed to rebuild wire while replacing chain with midline";
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
                    errorMessage = "ShapeFix_Wire returned null wire after midline replacement";
                    return false;
                }

                fixedWire.Orientation(oldWire.Orientation());
                wireToAdd = fixedWire;
            }

            builder.Add(rebuilt, wireToAdd);
        }

        if (!replacedChain) {
            errorMessage = "face does not contain a contiguous source chain";
            return false;
        }

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
            errorMessage = "ShapeFix_Face returned null face after midline replacement";
            return false;
        }

        rebuiltFace.Orientation(faceToRebuild.Orientation());
    }
    catch (Standard_Failure& failure) {
        errorMessage = "face rebuild with midline failed: " + failureText(failure);
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
    const std::vector<BRepCheckFailureInfo> existingFailures =
        collectBRepCheckFailures(beforeShape, &existingError);
    if (!existingError.empty()) {
        errorMessage = "existing BRepCheck failure collection failed: " + existingError;
        return false;
    }

    std::string repairedError;
    const std::vector<BRepCheckFailureInfo> repairedFailures =
        collectBRepCheckFailures(afterShape, &repairedError);
    if (!repairedError.empty()) {
        errorMessage = "repaired BRepCheck failure collection failed: " + repairedError;
        return false;
    }

    std::vector<BRepCheckFailureInfo> newFailures;
    for (const BRepCheckFailureInfo& failure : repairedFailures) {
        if (!isExistingBRepCheckFailure(failure, existingFailures)) {
            newFailures.push_back(failure);
        }
    }

    if (!newFailures.empty() && params_.rejectInvalidResult) {
        errorMessage = "new BRepCheck failures after closed narrow gap repair:"
            + describeBRepCheckFailures(newFailures);
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

bool HealClosedNarrowGap::findContiguousChainInWire(
    const std::vector<TopoDS_Edge>& wireEdges,
    const EdgeChain& chain,
    int& startIndex,
    bool& reversedOrder) const
{
    startIndex = -1;
    reversedOrder = false;
    if (wireEdges.empty() || chain.edges.empty() || chain.edges.size() > wireEdges.size()) {
        return false;
    }

    const int wireEdgeCount = static_cast<int>(wireEdges.size());
    const int chainEdgeCount = static_cast<int>(chain.edges.size());
    for (int start = 0; start < wireEdgeCount; ++start) {
        bool forwardMatch = true;
        for (int offset = 0; offset < chainEdgeCount; ++offset) {
            const int wireIndex = (start + offset) % wireEdgeCount;
            if (!wireEdges[static_cast<std::size_t>(wireIndex)].IsSame(
                    chain.edges[static_cast<std::size_t>(offset)])) {
                forwardMatch = false;
                break;
            }
        }

        if (forwardMatch) {
            startIndex = start;
            reversedOrder = false;
            return true;
        }

        bool reverseMatch = true;
        for (int offset = 0; offset < chainEdgeCount; ++offset) {
            const int chainIndex = chainEdgeCount - 1 - offset;
            const int wireIndex = (start + offset) % wireEdgeCount;
            if (!wireEdges[static_cast<std::size_t>(wireIndex)].IsSame(
                    chain.edges[static_cast<std::size_t>(chainIndex)])) {
                reverseMatch = false;
                break;
            }
        }

        if (reverseMatch) {
            startIndex = start;
            reversedOrder = true;
            return true;
        }
    }

    return false;
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
