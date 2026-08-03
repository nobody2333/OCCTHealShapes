#include "HealClosedNarrowGap.h"

#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
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
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_POnCurv.hxx>
#include <GCPnts_AbscissaPoint.hxx>
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
constexpr double kStationMergeTolerance =
    128.0 * std::numeric_limits<double>::epsilon();

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

bool isFinitePoint(const gp_Pnt& point)
{
    return std::isfinite(point.X())
        && std::isfinite(point.Y())
        && std::isfinite(point.Z());
}

std::string formatIndexList(const std::vector<std::size_t>& indices)
{
    std::ostringstream os;
    os << "[";
    for (std::size_t index = 0; index < indices.size(); ++index) {
        if (index > 0) {
            os << ",";
        }
        os << indices[index];
    }
    os << "]";
    return os.str();
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
    params_.widthSampleCount = std::max(2, params_.widthSampleCount);
    params_.minWidthSampleCount = std::max(2, params_.minWidthSampleCount);
    params_.maxWidthSampleCount = std::max(
        params_.minWidthSampleCount,
        params_.maxWidthSampleCount);
    if (!std::isfinite(params_.targetWidthSampleSpacing) ||
        params_.targetWidthSampleSpacing <= Precision::Confusion()) {
        params_.useAdaptiveWidthSampleCount = false;
    }
    params_.sampleCount = std::max(2, params_.sampleCount);
    params_.minSampleCount = std::max(2, params_.minSampleCount);
    params_.maxSampleCount = std::max(params_.minSampleCount, params_.maxSampleCount);
    if (params_.targetSampleSpacing <= Precision::Confusion()) {
        params_.useAdaptiveSampleCount = false;
    }
    params_.midlineRepairTolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
    params_.projectionVertexSnapTolerance = std::max(
        params_.projectionVertexSnapTolerance,
        Precision::Confusion());
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
           << ", widthSampleCount=" << params_.widthSampleCount
           << ", useAdaptiveWidthSampleCount="
           << (params_.useAdaptiveWidthSampleCount ? "true" : "false")
           << ", minWidthSampleCount=" << params_.minWidthSampleCount
           << ", maxWidthSampleCount=" << params_.maxWidthSampleCount
           << ", targetWidthSampleSpacing="
           << params_.targetWidthSampleSpacing
           << ", sampleCount=" << params_.sampleCount
           << ", useAdaptiveSampleCount="
           << (params_.useAdaptiveSampleCount ? "true" : "false")
           << ", minSampleCount=" << params_.minSampleCount
           << ", maxSampleCount=" << params_.maxSampleCount
           << ", targetSampleSpacing=" << params_.targetSampleSpacing
           << ", useMidlineCollapse=" << (params_.useMidlineCollapse ? "true" : "false")
           << ", midlineConstructionMode="
           << midlineConstructionModeName(params_.midlineConstructionMode)
           << ", projectionVertexSnapTolerance="
           << params_.projectionVertexSnapTolerance
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
                   << "\n  maxWidthAToB=" << candidate.maxWidthAToB
                   << "\n  maxWidthBToA=" << candidate.maxWidthBToA
                   << "\n  meanWidth=" << candidate.meanWidth
                   << "\n  widthSampleCount="
                   << candidate.usedWidthSampleCount
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
               << ", maxWidthAToB=" << candidate.maxWidthAToB
               << ", maxWidthBToA=" << candidate.maxWidthBToA
               << ", meanWidth=" << candidate.meanWidth
               << ", widthSampleCount="
               << candidate.usedWidthSampleCount
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

                std::string gapError;
                if (!evaluateBidirectionalGap(
                        candidate.chainA,
                        candidate.chainB,
                        candidate.maxWidthAToB,
                        candidate.maxWidthBToA,
                        candidate.meanWidth,
                        candidate.usedWidthSampleCount,
                        gapError)) {
                    continue;
                }

                candidate.maxWidth = std::max(
                    candidate.maxWidthAToB,
                    candidate.maxWidthBToA);
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
        !std::isfinite(chain.length) ||
        chain.length <= Precision::Confusion()) {
        return false;
    }

    try {
        const double targetLength = chain.length * clamp01(normalized);
        double accumulated = 0.0;
        for (std::size_t i = 0; i < chain.edges.size(); ++i) {
            const TopoDS_Edge& edge = chain.edges[i];
            const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
            if (!std::isfinite(edgeLength) ||
                edgeLength <= Precision::Confusion()) {
                return false;
            }

            const bool isLastEdge = i + 1 == chain.edges.size();
            if (!isLastEdge &&
                accumulated + edgeLength < targetLength) {
                accumulated += edgeLength;
                continue;
            }

            const double localLength = std::max(
                0.0,
                std::min(edgeLength, targetLength - accumulated));
            const double startParameter =
                BRep_Tool::Parameter(chain.vertices[i], edge);
            const double endParameter =
                BRep_Tool::Parameter(chain.vertices[i + 1], edge);
            if (!std::isfinite(startParameter) ||
                !std::isfinite(endParameter) ||
                std::abs(endParameter - startParameter) <=
                    Precision::PConfusion()) {
                return false;
            }

            // BRepAdaptor_Curve 会应用 edge 的 TopLoc_Location，投影点与源点位于同一全局坐标系。
            BRepAdaptor_Curve curve(edge);
            double parameter = startParameter;
            if (localLength >= edgeLength - Precision::Confusion()) {
                parameter = endParameter;
            }
            else if (localLength > Precision::Confusion()) {
                const double initialParameter =
                    startParameter
                    + (endParameter - startParameter)
                        * (localLength / edgeLength);
                const double signedAbscissa =
                    startParameter <= endParameter
                        ? localLength
                        : -localLength;
                GCPnts_AbscissaPoint abscissaPoint(
                    curve,
                    signedAbscissa,
                    startParameter,
                    initialParameter,
                    Precision::Confusion());
                if (!abscissaPoint.IsDone()) {
                    return false;
                }
                parameter = abscissaPoint.Parameter();
            }

            point = curve.Value(parameter);
            return isFinitePoint(point);
        }
    }
    catch (const Standard_Failure&) {
        return false;
    }

    return false;
}

bool HealClosedNarrowGap::findNearestPointOnChain(
    const gp_Pnt& sourcePoint,
    const EdgeChain& targetChain,
    std::size_t& targetEdgeIndex,
    double& targetParameter,
    gp_Pnt& targetPoint,
    double& distance,
    std::string& errorMessage) const
{
    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    targetEdgeIndex = invalidIndex;
    targetParameter = 0.0;
    targetPoint = gp_Pnt();
    distance = std::numeric_limits<double>::infinity();

    if (!isFinitePoint(sourcePoint) ||
        targetChain.edges.empty() ||
        targetChain.vertices.size() != targetChain.edges.size() + 1) {
        errorMessage = "nearest chain projection failure: target chain data is invalid";
        return false;
    }

    try {
        for (std::size_t edgeIndex = 0;
             edgeIndex < targetChain.edges.size();
             ++edgeIndex) {
            const TopoDS_Edge& edge = targetChain.edges[edgeIndex];
            if (edge.IsNull() ||
                targetChain.vertices[edgeIndex].IsNull() ||
                targetChain.vertices[edgeIndex + 1].IsNull()) {
                errorMessage = "nearest chain projection failure: target chain contains null topology";
                return false;
            }

            BRepAdaptor_Curve curve(edge);
            const double startParameter =
                BRep_Tool::Parameter(targetChain.vertices[edgeIndex], edge);
            const double endParameter =
                BRep_Tool::Parameter(targetChain.vertices[edgeIndex + 1], edge);
            const double lowerParameter = std::max(
                curve.FirstParameter(),
                std::min(startParameter, endParameter));
            const double upperParameter = std::min(
                curve.LastParameter(),
                std::max(startParameter, endParameter));
            if (!std::isfinite(startParameter) ||
                !std::isfinite(endParameter) ||
                !std::isfinite(lowerParameter) ||
                !std::isfinite(upperParameter) ||
                upperParameter - lowerParameter <= Precision::PConfusion()) {
                errorMessage = "nearest chain projection failure: target edge has an invalid trimmed range";
                return false;
            }

            auto considerProjection =
                [&](double parameter, const gp_Pnt& candidatePoint) {
                    if (!std::isfinite(parameter) ||
                        !isFinitePoint(candidatePoint) ||
                        parameter < lowerParameter - Precision::PConfusion() ||
                        parameter > upperParameter + Precision::PConfusion()) {
                        return;
                    }

                    const double candidateDistance =
                        sourcePoint.Distance(candidatePoint);
                    if (!std::isfinite(candidateDistance)) {
                        return;
                    }
                    if (targetEdgeIndex == invalidIndex ||
                        candidateDistance + Precision::Confusion() < distance) {
                        distance = candidateDistance;
                        targetEdgeIndex = edgeIndex;
                        targetParameter = parameter;
                        targetPoint = candidatePoint;
                    }
                };

            // 极值求解限制在实际裁剪参数范围内，并显式加入两个端点作为兜底候选。
            considerProjection(lowerParameter, curve.Value(lowerParameter));
            considerProjection(upperParameter, curve.Value(upperParameter));
            Extrema_ExtPC extrema(
                sourcePoint,
                curve,
                lowerParameter,
                upperParameter,
                Precision::Confusion());
            if (extrema.IsDone()) {
                for (int solutionIndex = 1;
                     solutionIndex <= extrema.NbExt();
                     ++solutionIndex) {
                    const Extrema_POnCurv& pointOnCurve =
                        extrema.Point(solutionIndex);
                    considerProjection(
                        pointOnCurve.Parameter(),
                        pointOnCurve.Value());
                }
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "nearest chain projection failure: "
            + failureText(failure);
        return false;
    }

    if (targetEdgeIndex == invalidIndex ||
        !std::isfinite(targetParameter) ||
        !isFinitePoint(targetPoint) ||
        !std::isfinite(distance)) {
        errorMessage = "nearest chain projection failure: no valid projection on trimmed target edges";
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::evaluateBidirectionalGap(
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    double& maxWidthAToB,
    double& maxWidthBToA,
    double& meanWidth,
    int& usedSampleCount,
    std::string& errorMessage) const
{
    maxWidthAToB = 0.0;
    maxWidthBToA = 0.0;
    meanWidth = 0.0;
    usedSampleCount = computeWidthSampleCount(
        chainA.length,
        chainB.length);
    const int sampleCount = usedSampleCount;
    double distanceSum = 0.0;
    int distanceCount = 0;

    auto accumulateDirectedDistances =
        [&](const EdgeChain& sourceChain,
            const EdgeChain& targetChain,
            double& directedMaximum,
            const char* directionName) {
            for (int sampleIndex = 0;
                 sampleIndex < sampleCount;
                 ++sampleIndex) {
                const double normalized =
                    static_cast<double>(sampleIndex)
                    / static_cast<double>(sampleCount - 1);
                gp_Pnt sourcePoint;
                if (!sampleChainPoint(
                        sourceChain,
                        normalized,
                        sourcePoint)) {
                    errorMessage = std::string(
                        "bidirectional gap evaluation failure: failed to sample ")
                        + directionName;
                    return false;
                }

                std::size_t targetEdgeIndex = 0;
                double targetParameter = 0.0;
                gp_Pnt targetPoint;
                double nearestDistance = 0.0;
                std::string projectionError;
                if (!findNearestPointOnChain(
                        sourcePoint,
                        targetChain,
                        targetEdgeIndex,
                        targetParameter,
                        targetPoint,
                        nearestDistance,
                        projectionError)) {
                    errorMessage = "bidirectional gap evaluation failure: "
                        + projectionError;
                    return false;
                }

                directedMaximum = std::max(
                    directedMaximum,
                    nearestDistance);
                distanceSum += nearestDistance;
                ++distanceCount;
            }
            return true;
        };

    // 两个方向都必须覆盖，避免只检查 A 靠近 B，却遗漏 B 的局部外张。
    if (!accumulateDirectedDistances(
            chainA,
            chainB,
            maxWidthAToB,
            "chainA") ||
        !accumulateDirectedDistances(
            chainB,
            chainA,
            maxWidthBToA,
            "chainB")) {
        return false;
    }

    if (distanceCount != 2 * sampleCount ||
        !std::isfinite(distanceSum) ||
        !std::isfinite(maxWidthAToB) ||
        !std::isfinite(maxWidthBToA)) {
        errorMessage = "bidirectional gap evaluation failure: distance statistics are invalid";
        return false;
    }

    meanWidth = distanceSum / static_cast<double>(distanceCount);
    if (!std::isfinite(meanWidth)) {
        errorMessage = "bidirectional gap evaluation failure: mean distance is invalid";
        return false;
    }
    return true;
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
    if (candidate.chainA.vertices.size() != candidate.chainA.edges.size() + 1 ||
        candidate.chainB.vertices.size() != candidate.chainB.edges.size() + 1) {
        errorMessage = "candidate chain vertex count is invalid";
        return false;
    }
    if (candidate.chainA.ownerFaces.size() != candidate.chainA.edges.size() ||
        candidate.chainB.ownerFaces.size() != candidate.chainB.edges.size()) {
        errorMessage = "candidate chain owner-face count is invalid";
        return false;
    }
    if (!std::isfinite(candidate.chainA.length) ||
        !std::isfinite(candidate.chainB.length) ||
        candidate.chainA.length <= Precision::Confusion() ||
        candidate.chainB.length <= Precision::Confusion()) {
        errorMessage = "candidate chain total length is invalid";
        return false;
    }
    if (!candidate.chainA.vertices.front().IsSame(candidate.chainB.vertices.front()) ||
        !candidate.chainA.vertices.back().IsSame(candidate.chainB.vertices.back())) {
        errorMessage = "candidate chains do not share the same endpoints";
        return false;
    }
    for (const TopoDS_Face& ownerFace : candidate.chainA.ownerFaces) {
        if (ownerFace.IsNull()) {
            errorMessage = "chainA contains an edge without owner face";
            return false;
        }
    }
    for (const TopoDS_Face& ownerFace : candidate.chainB.ownerFaces) {
        if (ownerFace.IsNull()) {
            errorMessage = "chainB contains an edge without owner face";
            return false;
        }
    }

    MidlineConstructionMode constructionMode =
        MidlineConstructionMode::Auto;
    if (!resolveMidlineConstructionMode(
            candidate,
            constructionMode,
            errorMessage)) {
        return false;
    }

    if (constructionMode == MidlineConstructionMode::LegacyOneToOne) {
        return buildMidlineChainOneToOne(candidate, midlineData, errorMessage);
    }
    if (constructionMode ==
        MidlineConstructionMode::ProjectionCorrespondence) {
        return buildMidlineChainByProjectionCorrespondence(
            candidate,
            midlineData,
            errorMessage);
    }

    errorMessage = "resolved midline construction mode is invalid";
    return false;
}

bool HealClosedNarrowGap::canUseLegacyOneToOnePath(
    const ClosedGapCandidate& candidate) const
{
    if (candidate.chainA.edges.size() != candidate.chainB.edges.size() ||
        candidate.chainA.vertices.size() != candidate.chainB.vertices.size() ||
        candidate.chainA.vertices.size() != candidate.chainA.edges.size() + 1) {
        return false;
    }

    try {
        for (std::size_t vertexIndex = 1;
             vertexIndex + 1 < candidate.chainA.vertices.size();
             ++vertexIndex) {
            const TopoDS_Vertex& vertexA = candidate.chainA.vertices[vertexIndex];
            const TopoDS_Vertex& vertexB = candidate.chainB.vertices[vertexIndex];
            if (vertexA.IsNull() || vertexB.IsNull()) {
                return false;
            }

            const gp_Pnt pointA = BRep_Tool::Pnt(vertexA);
            const gp_Pnt pointB = BRep_Tool::Pnt(vertexB);
            if (!isFinitePoint(pointA) || !isFinitePoint(pointB) ||
                pointA.Distance(pointB) >
                    params_.maxGapWidth + params_.midlineRepairTolerance) {
                return false;
            }
        }
    }
    catch (Standard_Failure&) {
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::resolveMidlineConstructionMode(
    const ClosedGapCandidate& candidate,
    MidlineConstructionMode& mode,
    std::string& errorMessage) const
{
    switch (params_.midlineConstructionMode) {
    case MidlineConstructionMode::Auto:
        // 自动模式保留原有行为：可安全一一配对时继续使用旧路径。
        mode = canUseLegacyOneToOnePath(candidate)
            ? MidlineConstructionMode::LegacyOneToOne
            : MidlineConstructionMode::ProjectionCorrespondence;
        return true;

    case MidlineConstructionMode::LegacyOneToOne:
        if (!canUseLegacyOneToOnePath(candidate)) {
            errorMessage =
                "forced LegacyOneToOne mode requires equal source edge counts "
                "and safely paired internal vertices";
            return false;
        }
        mode = MidlineConstructionMode::LegacyOneToOne;
        return true;

    case MidlineConstructionMode::ProjectionCorrespondence:
        // 投影对应是通用路径，因此也允许显式用于对称候选。
        mode = MidlineConstructionMode::ProjectionCorrespondence;
        return true;

    default:
        errorMessage = "unknown midline construction mode";
        return false;
    }
}

bool HealClosedNarrowGap::buildMidlineChainOneToOne(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    midlineData.constructionMode = MidlineConstructionMode::LegacyOneToOne;

    try {
        const std::size_t edgeCount = candidate.chainA.edges.size();
        std::vector<TopoDS_Vertex> midVertices;
        midVertices.reserve(edgeCount + 1);

        BRep_Builder builder;
        for (std::size_t vertexIndex = 0; vertexIndex <= edgeCount; ++vertexIndex) {
            if (vertexIndex == 0) {
                // 两侧起始尖点已经拓扑重合，直接沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.front());
                continue;
            }
            if (vertexIndex == edgeCount) {
                // 两侧结束尖点同样沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.back());
                continue;
            }

            const TopoDS_Vertex& vertexA = candidate.chainA.vertices[vertexIndex];
            const TopoDS_Vertex& vertexB = candidate.chainB.vertices[vertexIndex];
            if (vertexA.IsNull() || vertexB.IsNull()) {
                errorMessage = "target vertex construction failure: "
                    "candidate contains a null corresponding internal vertex";
                return false;
            }

            const gp_Pnt pointA = BRep_Tool::Pnt(vertexA);
            const gp_Pnt pointB = BRep_Tool::Pnt(vertexB);
            if (!isFinitePoint(pointA) || !isFinitePoint(pointB)) {
                errorMessage = "target vertex construction failure: "
                    "candidate contains a non-finite internal vertex";
                return false;
            }
            const double vertexGap = pointA.Distance(pointB);
            if (vertexGap > params_.maxGapWidth + params_.midlineRepairTolerance) {
                errorMessage = "target vertex construction failure: "
                    "corresponding internal chain vertices are too far apart";
                return false;
            }

            const gp_Pnt midpoint(
                0.5 * (pointA.X() + pointB.X()),
                0.5 * (pointA.Y() + pointB.Y()),
                0.5 * (pointA.Z() + pointB.Z()));

            // 新顶点还会替换两侧旧顶点，因此容差必须覆盖旧顶点到中点的移动距离。
            const double vertexTolerance = std::max({
                Precision::Confusion(),
                params_.midlineRepairTolerance,
                BRep_Tool::Tolerance(vertexA),
                BRep_Tool::Tolerance(vertexB),
                midpoint.Distance(pointA) + Precision::Confusion(),
                midpoint.Distance(pointB) + Precision::Confusion()
            });

            TopoDS_Vertex midVertex;
            builder.MakeVertex(midVertex, midpoint, vertexTolerance);
            if (midVertex.IsNull()) {
                errorMessage = "target vertex construction failure: "
                    "failed to create midline vertex";
                return false;
            }
            midVertices.push_back(midVertex);
        }

        midlineData.chain.vertices = midVertices;
        midlineData.replacementsForChainA.resize(edgeCount);
        midlineData.replacementsForChainB.resize(edgeCount);
        midlineData.targetVertexIndicesForChainA.resize(edgeCount + 1);
        midlineData.targetVertexIndicesForChainB.resize(edgeCount + 1);
        for (std::size_t vertexIndex = 0; vertexIndex <= edgeCount; ++vertexIndex) {
            midlineData.targetVertexIndicesForChainA[vertexIndex] = vertexIndex;
            midlineData.targetVertexIndicesForChainB[vertexIndex] = vertexIndex;
        }

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
                if (!isFinitePoint(sampledPoints.back())) {
                    errorMessage = std::string(chainName)
                        + " edge produced a non-finite sample point";
                    return false;
                }
            }
            return true;
        };

        for (std::size_t edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
            const double lengthA = GeometryUtils::ComputeEdgeLength(candidate.chainA.edges[edgeIndex]);
            const double lengthB = GeometryUtils::ComputeEdgeLength(candidate.chainB.edges[edgeIndex]);
            if (!std::isfinite(lengthA) || !std::isfinite(lengthB) ||
                lengthA <= Precision::Confusion() ||
                lengthB <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: "
                    "candidate contains an invalid source edge pair";
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
                errorMessage = "target edge construction failure: "
                    "midline curve guidance points collapse at the segment endpoint";
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
                errorMessage = "target edge construction failure: "
                    "failed to interpolate a midline curve from sampled midpoints";
                return false;
            }

            Handle(Geom_Curve) midlineCurve = interpolator.Curve();
            BRepBuilderAPI_MakeEdge edgeMaker(
                midlineCurve,
                midVertices[edgeIndex],
                midVertices[edgeIndex + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "target edge construction failure: "
                    "failed to create midline edge";
                return false;
            }

            const TopoDS_Edge edge = edgeMaker.Edge();
            if (edge.IsNull()) {
                errorMessage = "target edge construction failure: "
                    "created midline edge is null";
                return false;
            }
            const double midlineEdgeLength = GeometryUtils::ComputeEdgeLength(edge);
            if (!std::isfinite(midlineEdgeLength) ||
                midlineEdgeLength <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: "
                    "created midline edge has invalid length";
                return false;
            }
            midlineData.chain.edges.push_back(edge);
            midlineData.chain.length += midlineEdgeLength;
            const std::size_t midlineEdgeIndex = midlineData.chain.edges.size() - 1;

            // source edge 与中线 edge 一一对应，采样密度不再改变拓扑分段数量。
            midlineData.replacementsForChainA[edgeIndex].push_back(midlineEdgeIndex);
            midlineData.replacementsForChainB[edgeIndex].push_back(midlineEdgeIndex);
        }

        if (midlineData.chain.edges.empty() ||
            !std::isfinite(midlineData.chain.length) ||
            midlineData.chain.length <= Precision::Confusion()) {
            errorMessage = "midline chain is empty";
            return false;
        }

        if (!validateReplacementMapping(
                midlineData.replacementsForChainA,
                candidate.chainA.edges.size(),
                midlineData.chain.edges.size(),
                "chainA",
                errorMessage) ||
            !validateReplacementMapping(
                midlineData.replacementsForChainB,
                candidate.chainB.edges.size(),
                midlineData.chain.edges.size(),
                "chainB",
                errorMessage)) {
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "legacy one-to-one midline construction failed: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::buildNormalizedVertexStations(
    const EdgeChain& chain,
    const char* chainName,
    ChainStationData& stationData,
    std::string& errorMessage) const
{
    stationData = ChainStationData();
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (chain.edges.empty() ||
        chain.vertices.size() != chain.edges.size() + 1 ||
        chain.ownerFaces.size() != chain.edges.size()) {
        errorMessage = "station construction failure: " + name
            + " topology counts are inconsistent";
        return false;
    }
    if (!std::isfinite(chain.length) ||
        chain.length <= Precision::Confusion()) {
        errorMessage = "station construction failure: " + name
            + " total length is invalid";
        return false;
    }

    stationData.edgeLengths.reserve(chain.edges.size());
    for (const TopoDS_Edge& edge : chain.edges) {
        const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
        if (!std::isfinite(edgeLength) ||
            edgeLength <= Precision::Confusion()) {
            errorMessage = "station construction failure: " + name
                + " contains an invalid source edge length";
            return false;
        }

        stationData.edgeLengths.push_back(edgeLength);
        stationData.totalLength += edgeLength;
        if (!std::isfinite(stationData.totalLength)) {
            errorMessage = "station construction failure: " + name
                + " accumulated length is not finite";
            return false;
        }
    }

    if (stationData.totalLength <= Precision::Confusion()) {
        errorMessage = "station construction failure: " + name
            + " computed total length is invalid";
        return false;
    }

    stationData.vertexStations.reserve(chain.vertices.size());
    stationData.vertexStations.push_back(0.0);
    double accumulatedLength = 0.0;
    for (std::size_t edgeIndex = 0;
         edgeIndex + 1 < stationData.edgeLengths.size();
         ++edgeIndex) {
        accumulatedLength += stationData.edgeLengths[edgeIndex];
        const double station = accumulatedLength / stationData.totalLength;
        if (!std::isfinite(station) ||
            station <= stationData.vertexStations.back() ||
            station >= 1.0) {
            errorMessage = "station construction failure: " + name
                + " vertex stations are not strictly increasing";
            return false;
        }
        stationData.vertexStations.push_back(station);
    }
    stationData.vertexStations.push_back(1.0);

    if (stationData.vertexStations.size() != chain.vertices.size()) {
        errorMessage = "station construction failure: " + name
            + " vertex station count is inconsistent";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::projectPointToChain(
    const gp_Pnt& sourcePoint,
    const EdgeChain& targetChain,
    const ChainStationData& targetStationData,
    ChainProjectionResult& projection,
    std::string& errorMessage) const
{
    projection = ChainProjectionResult();
    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    if (!isFinitePoint(sourcePoint) ||
        targetChain.edges.empty() ||
        targetChain.vertices.size() != targetChain.edges.size() + 1 ||
        targetStationData.edgeLengths.size() != targetChain.edges.size() ||
        targetStationData.vertexStations.size() != targetChain.vertices.size() ||
        !std::isfinite(targetStationData.totalLength) ||
        targetStationData.totalLength <= Precision::Confusion()) {
        errorMessage = "chain projection failure: target chain data is invalid";
        return false;
    }

    std::size_t bestEdgeIndex = invalidIndex;
    double bestParameter = 0.0;
    gp_Pnt bestPoint;
    double bestDistance = std::numeric_limits<double>::infinity();
    std::string nearestError;
    if (!findNearestPointOnChain(
            sourcePoint,
            targetChain,
            bestEdgeIndex,
            bestParameter,
            bestPoint,
            bestDistance,
            nearestError)) {
        errorMessage = "chain projection failure: " + nearestError;
        return false;
    }

    try {
        const double allowedDistance =
            params_.maxGapWidth + params_.midlineRepairTolerance;
        if (bestDistance > allowedDistance) {
            errorMessage = "chain projection failure: nearest point exceeds allowed gap width";
            return false;
        }

        const TopoDS_Edge& bestEdge =
            targetChain.edges[bestEdgeIndex];
        BRepAdaptor_Curve bestCurve(bestEdge);
        const double edgeStartParameter = BRep_Tool::Parameter(
            targetChain.vertices[bestEdgeIndex],
            bestEdge);
        double localLength = GCPnts_AbscissaPoint::Length(
            bestCurve,
            edgeStartParameter,
            bestParameter,
            Precision::Confusion());
        if (!std::isfinite(localLength) ||
            localLength < -Precision::Confusion() ||
            localLength >
                targetStationData.edgeLengths[bestEdgeIndex]
                    + Precision::Confusion()) {
            errorMessage = "chain projection failure: projected arc length is invalid";
            return false;
        }
        localLength = std::max(
            0.0,
            std::min(
                targetStationData.edgeLengths[bestEdgeIndex],
                localLength));

        double accumulatedLength = 0.0;
        for (std::size_t edgeIndex = 0;
             edgeIndex < bestEdgeIndex;
             ++edgeIndex) {
            accumulatedLength += targetStationData.edgeLengths[edgeIndex];
        }
        const double targetStation = clamp01(
            (accumulatedLength + localLength)
            / targetStationData.totalLength);
        if (!std::isfinite(targetStation)) {
            errorMessage = "chain projection failure: projected chain station is invalid";
            return false;
        }

        projection.projectedPoint = bestPoint;
        projection.rawProjectedPoint = bestPoint;
        projection.targetEdgeIndex = bestEdgeIndex;
        projection.targetParameter = bestParameter;
        projection.targetStation = targetStation;
        projection.rawTargetStation = targetStation;
        projection.distance = bestDistance;
        projection.rawDistance = bestDistance;

        const double snapTolerance = std::min(
            params_.projectionVertexSnapTolerance,
            params_.maxGapWidth);
        const double stationSnapTolerance = std::min(
            0.25,
            snapTolerance / std::max(
                targetStationData.totalLength,
                Precision::Confusion()));

        std::size_t snappedVertexIndex = invalidIndex;
        double snappedVertexDistance =
            std::numeric_limits<double>::infinity();
        const std::size_t endpointIndices[2] = {
            bestEdgeIndex,
            bestEdgeIndex + 1
        };
        for (const std::size_t vertexIndex : endpointIndices) {
            const gp_Pnt vertexPoint =
                BRep_Tool::Pnt(targetChain.vertices[vertexIndex]);
            const double geometricDistance =
                bestPoint.Distance(vertexPoint);
            const double stationDistance = std::abs(
                targetStation
                - targetStationData.vertexStations[vertexIndex]);
            if (!std::isfinite(geometricDistance) ||
                geometricDistance > snapTolerance ||
                stationDistance > stationSnapTolerance) {
                continue;
            }

            if (snappedVertexIndex == invalidIndex ||
                geometricDistance + Precision::Confusion()
                    < snappedVertexDistance ||
                (std::abs(geometricDistance - snappedVertexDistance)
                        <= Precision::Confusion() &&
                    vertexIndex < snappedVertexIndex)) {
                snappedVertexIndex = vertexIndex;
                snappedVertexDistance = geometricDistance;
            }
        }

        // 吸附容差独立于细缝宽度，并同时检查链站点，避免吸到空间近但拓扑远的顶点。
        if (snappedVertexIndex != invalidIndex) {
            projection.projectedPoint =
                BRep_Tool::Pnt(targetChain.vertices[snappedVertexIndex]);
            projection.targetStation =
                targetStationData.vertexStations[snappedVertexIndex];
            projection.distance =
                sourcePoint.Distance(projection.projectedPoint);
            projection.snappedTargetVertexIndex =
                snappedVertexIndex;
            projection.snappedToVertex = true;
            if (!std::isfinite(projection.distance) ||
                projection.distance > allowedDistance) {
                errorMessage = "chain projection failure: snapped vertex exceeds allowed gap width";
                return false;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "chain projection failure: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::selectSnappedVertexPairs(
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    const std::vector<ChainProjectionResult>& projectionsFromA,
    const std::vector<ChainProjectionResult>& projectionsFromB,
    std::vector<SnappedVertexPair>& selectedPairs,
    std::string& errorMessage) const
{
    selectedPairs.clear();
    if (projectionsFromA.size() != chainA.vertices.size() ||
        projectionsFromB.size() != chainB.vertices.size()) {
        errorMessage = "projection conflict failure: projection counts are inconsistent";
        return false;
    }

    struct PairProposal {
        std::size_t vertexIndexA = 0;
        std::size_t vertexIndexB = 0;
        double distance = 0.0;
    };

    std::vector<PairProposal> proposals;
    auto addProposal = [&](std::size_t vertexIndexA,
                           std::size_t vertexIndexB,
                           double distance) {
        if (vertexIndexA == 0 ||
            vertexIndexA + 1 >= chainA.vertices.size() ||
            vertexIndexB == 0 ||
            vertexIndexB + 1 >= chainB.vertices.size() ||
            !std::isfinite(distance)) {
            return;
        }

        for (PairProposal& proposal : proposals) {
            if (proposal.vertexIndexA == vertexIndexA &&
                proposal.vertexIndexB == vertexIndexB) {
                proposal.distance = std::min(proposal.distance, distance);
                return;
            }
        }

        PairProposal proposal;
        proposal.vertexIndexA = vertexIndexA;
        proposal.vertexIndexB = vertexIndexB;
        proposal.distance = distance;
        proposals.push_back(proposal);
    };

    for (std::size_t vertexIndexA = 1;
         vertexIndexA + 1 < chainA.vertices.size();
         ++vertexIndexA) {
        const ChainProjectionResult& projection =
            projectionsFromA[vertexIndexA];
        if (projection.snappedToVertex) {
            addProposal(
                vertexIndexA,
                projection.snappedTargetVertexIndex,
                projection.distance);
        }
    }
    for (std::size_t vertexIndexB = 1;
         vertexIndexB + 1 < chainB.vertices.size();
         ++vertexIndexB) {
        const ChainProjectionResult& projection =
            projectionsFromB[vertexIndexB];
        if (projection.snappedToVertex) {
            addProposal(
                projection.snappedTargetVertexIndex,
                vertexIndexB,
                projection.distance);
        }
    }

    std::sort(
        proposals.begin(),
        proposals.end(),
        [](const PairProposal& lhs, const PairProposal& rhs) {
            if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
            }
            if (lhs.vertexIndexA != rhs.vertexIndexA) {
                return lhs.vertexIndexA < rhs.vertexIndexA;
            }
            return lhs.vertexIndexB < rhs.vertexIndexB;
        });

    std::vector<bool> pairedVerticesA(chainA.vertices.size(), false);
    std::vector<bool> pairedVerticesB(chainB.vertices.size(), false);
    for (const PairProposal& proposal : proposals) {
        if (pairedVerticesA[proposal.vertexIndexA] ||
            pairedVerticesB[proposal.vertexIndexB]) {
            continue;
        }

        // 一个已有顶点只能配对一次；按距离和顶点下标排序保证冲突处理确定性。
        pairedVerticesA[proposal.vertexIndexA] = true;
        pairedVerticesB[proposal.vertexIndexB] = true;
        SnappedVertexPair selectedPair;
        selectedPair.sourceVertexIndexA = proposal.vertexIndexA;
        selectedPair.sourceVertexIndexB = proposal.vertexIndexB;
        selectedPair.distance = proposal.distance;
        selectedPairs.push_back(selectedPair);
    }

    std::sort(
        selectedPairs.begin(),
        selectedPairs.end(),
        [](const SnappedVertexPair& lhs,
           const SnappedVertexPair& rhs) {
            if (lhs.sourceVertexIndexA != rhs.sourceVertexIndexA) {
                return lhs.sourceVertexIndexA < rhs.sourceVertexIndexA;
            }
            return lhs.sourceVertexIndexB < rhs.sourceVertexIndexB;
        });
    return true;
}

bool HealClosedNarrowGap::buildProjectionCorrespondenceNodes(
    const ClosedGapCandidate& candidate,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    std::vector<ProjectionCorrespondenceNode>& nodes,
    std::string& errorMessage) const
{
    nodes.clear();
    if (candidate.chainA.vertices.size() !=
            chainAStations.vertexStations.size() ||
        candidate.chainB.vertices.size() !=
            chainBStations.vertexStations.size()) {
        errorMessage = "correspondence construction failure: source station counts are inconsistent";
        return false;
    }

    try {
        for (const TopoDS_Vertex& vertex : candidate.chainA.vertices) {
            if (vertex.IsNull() || !isFinitePoint(BRep_Tool::Pnt(vertex))) {
                errorMessage = "correspondence construction failure: chainA contains an invalid vertex";
                return false;
            }
        }
        for (const TopoDS_Vertex& vertex : candidate.chainB.vertices) {
            if (vertex.IsNull() || !isFinitePoint(BRep_Tool::Pnt(vertex))) {
                errorMessage = "correspondence construction failure: chainB contains an invalid vertex";
                return false;
            }
        }

        std::vector<ChainProjectionResult> projectionsFromA(
            candidate.chainA.vertices.size());
        std::vector<ChainProjectionResult> projectionsFromB(
            candidate.chainB.vertices.size());

        // 必须双向投影，否则只存在于另一侧的内部拓扑顶点会被遗漏。
        for (std::size_t vertexIndexA = 1;
             vertexIndexA + 1 < candidate.chainA.vertices.size();
             ++vertexIndexA) {
            if (!projectPointToChain(
                    BRep_Tool::Pnt(
                        candidate.chainA.vertices[vertexIndexA]),
                    candidate.chainB,
                    chainBStations,
                    projectionsFromA[vertexIndexA],
                    errorMessage)) {
                errorMessage = "chainA vertex projection failed: " + errorMessage;
                return false;
            }
        }
        for (std::size_t vertexIndexB = 1;
             vertexIndexB + 1 < candidate.chainB.vertices.size();
             ++vertexIndexB) {
            if (!projectPointToChain(
                    BRep_Tool::Pnt(
                        candidate.chainB.vertices[vertexIndexB]),
                    candidate.chainA,
                    chainAStations,
                    projectionsFromB[vertexIndexB],
                    errorMessage)) {
                errorMessage = "chainB vertex projection failed: " + errorMessage;
                return false;
            }
        }

        std::vector<SnappedVertexPair> selectedPairs;
        if (!selectSnappedVertexPairs(
                candidate.chainA,
                candidate.chainB,
                projectionsFromA,
                projectionsFromB,
                selectedPairs,
                errorMessage)) {
            return false;
        }

        std::vector<bool> pairedVerticesA(
            candidate.chainA.vertices.size(),
            false);
        std::vector<bool> pairedVerticesB(
            candidate.chainB.vertices.size(),
            false);

        ProjectionCorrespondenceNode startNode;
        startNode.stationA = 0.0;
        startNode.stationB = 0.0;
        startNode.pointA =
            BRep_Tool::Pnt(candidate.chainA.vertices.front());
        startNode.pointB =
            BRep_Tool::Pnt(candidate.chainB.vertices.front());
        startNode.endpoint = true;
        nodes.push_back(startNode);

        for (const SnappedVertexPair& pair : selectedPairs) {
            if (pair.sourceVertexIndexA == 0 ||
                pair.sourceVertexIndexA + 1 >=
                    candidate.chainA.vertices.size() ||
                pair.sourceVertexIndexB == 0 ||
                pair.sourceVertexIndexB + 1 >=
                    candidate.chainB.vertices.size()) {
                errorMessage = "correspondence construction failure: snapped pair index is invalid";
                return false;
            }

            pairedVerticesA[pair.sourceVertexIndexA] = true;
            pairedVerticesB[pair.sourceVertexIndexB] = true;

            ProjectionCorrespondenceNode node;
            node.stationA =
                chainAStations.vertexStations[
                    pair.sourceVertexIndexA];
            node.stationB =
                chainBStations.vertexStations[
                    pair.sourceVertexIndexB];
            node.pointA = BRep_Tool::Pnt(
                candidate.chainA.vertices[
                    pair.sourceVertexIndexA]);
            node.pointB = BRep_Tool::Pnt(
                candidate.chainB.vertices[
                    pair.sourceVertexIndexB]);
            node.sourceVertexIndexA =
                pair.sourceVertexIndexA;
            node.sourceVertexIndexB =
                pair.sourceVertexIndexB;
            node.projectionDistance = pair.distance;
            node.hasSourceVertexA = true;
            node.hasSourceVertexB = true;
            node.snappedPair = true;
            nodes.push_back(node);
        }

        for (std::size_t vertexIndexA = 1;
             vertexIndexA + 1 < candidate.chainA.vertices.size();
             ++vertexIndexA) {
            if (pairedVerticesA[vertexIndexA]) {
                continue;
            }

            const ChainProjectionResult& projection =
                projectionsFromA[vertexIndexA];
            ProjectionCorrespondenceNode node;
            node.stationA =
                chainAStations.vertexStations[vertexIndexA];
            node.stationB = projection.rawTargetStation;
            node.pointA = BRep_Tool::Pnt(
                candidate.chainA.vertices[vertexIndexA]);
            node.pointB = projection.rawProjectedPoint;
            node.sourceVertexIndexA = vertexIndexA;
            node.projectionDistance = projection.rawDistance;
            node.hasSourceVertexA = true;
            nodes.push_back(node);
        }

        for (std::size_t vertexIndexB = 1;
             vertexIndexB + 1 < candidate.chainB.vertices.size();
             ++vertexIndexB) {
            if (pairedVerticesB[vertexIndexB]) {
                continue;
            }

            const ChainProjectionResult& projection =
                projectionsFromB[vertexIndexB];
            ProjectionCorrespondenceNode node;
            node.stationA = projection.rawTargetStation;
            node.stationB =
                chainBStations.vertexStations[vertexIndexB];
            node.pointA = projection.rawProjectedPoint;
            node.pointB = BRep_Tool::Pnt(
                candidate.chainB.vertices[vertexIndexB]);
            node.sourceVertexIndexB = vertexIndexB;
            node.projectionDistance = projection.rawDistance;
            node.hasSourceVertexB = true;
            nodes.push_back(node);
        }

        ProjectionCorrespondenceNode endNode;
        endNode.stationA = 1.0;
        endNode.stationB = 1.0;
        endNode.pointA =
            BRep_Tool::Pnt(candidate.chainA.vertices.back());
        endNode.pointB =
            BRep_Tool::Pnt(candidate.chainB.vertices.back());
        endNode.endpoint = true;
        nodes.push_back(endNode);
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "correspondence construction failure: "
            + failureText(failure);
        return false;
    }

    return mergeProjectionCorrespondenceNodes(
        nodes,
        candidate.chainA,
        candidate.chainB,
        chainAStations,
        chainBStations,
        errorMessage);
}

bool HealClosedNarrowGap::mergeProjectionCorrespondenceNodes(
    std::vector<ProjectionCorrespondenceNode>& nodes,
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    std::string& errorMessage) const
{
    if (nodes.size() < 2) {
        errorMessage = "correspondence merge failure: node sequence is incomplete";
        return false;
    }

    const double snapTolerance = std::min(
        params_.projectionVertexSnapTolerance,
        params_.maxGapWidth);
    const double stationMergeToleranceA = std::min(
        0.05,
        std::max(
            kStationMergeTolerance,
            snapTolerance / std::max(
                chainAStations.totalLength,
                Precision::Confusion())));
    const double stationMergeToleranceB = std::min(
        0.05,
        std::max(
            kStationMergeTolerance,
            snapTolerance / std::max(
                chainBStations.totalLength,
                Precision::Confusion())));

    auto midpointOfNode =
        [](const ProjectionCorrespondenceNode& node) {
            return gp_Pnt(
                0.5 * (node.pointA.X() + node.pointB.X()),
                0.5 * (node.pointA.Y() + node.pointB.Y()),
                0.5 * (node.pointA.Z() + node.pointB.Z()));
        };

    std::sort(
        nodes.begin(),
        nodes.end(),
        [](const ProjectionCorrespondenceNode& lhs,
           const ProjectionCorrespondenceNode& rhs) {
            if (lhs.stationA != rhs.stationA) {
                return lhs.stationA < rhs.stationA;
            }
            if (lhs.stationB != rhs.stationB) {
                return lhs.stationB < rhs.stationB;
            }
            if (lhs.sourceVertexIndexA != rhs.sourceVertexIndexA) {
                return lhs.sourceVertexIndexA < rhs.sourceVertexIndexA;
            }
            return lhs.sourceVertexIndexB < rhs.sourceVertexIndexB;
        });

    std::vector<ProjectionCorrespondenceNode> mergedNodes;
    for (const ProjectionCorrespondenceNode& incomingNode : nodes) {
        if (!std::isfinite(incomingNode.stationA) ||
            !std::isfinite(incomingNode.stationB) ||
            incomingNode.stationA < 0.0 ||
            incomingNode.stationA > 1.0 ||
            incomingNode.stationB < 0.0 ||
            incomingNode.stationB > 1.0 ||
            !isFinitePoint(incomingNode.pointA) ||
            !isFinitePoint(incomingNode.pointB)) {
            errorMessage = "correspondence merge failure: node contains invalid geometry";
            return false;
        }

        std::size_t mergeIndex =
            std::numeric_limits<std::size_t>::max();
        for (std::size_t existingIndex = mergedNodes.size();
             existingIndex > 0;
             --existingIndex) {
            const ProjectionCorrespondenceNode& existingNode =
                mergedNodes[existingIndex - 1];
            if (incomingNode.stationA - existingNode.stationA >
                stationMergeToleranceA) {
                break;
            }
            if (std::abs(incomingNode.stationA - existingNode.stationA) >
                    stationMergeToleranceA ||
                std::abs(incomingNode.stationB - existingNode.stationB) >
                    stationMergeToleranceB ||
                midpointOfNode(incomingNode).Distance(
                    midpointOfNode(existingNode)) > snapTolerance) {
                continue;
            }
            if (incomingNode.endpoint != existingNode.endpoint) {
                continue;
            }
            mergeIndex = existingIndex - 1;
            break;
        }

        if (mergeIndex ==
            std::numeric_limits<std::size_t>::max()) {
            mergedNodes.push_back(incomingNode);
            continue;
        }

        ProjectionCorrespondenceNode& existingNode =
            mergedNodes[mergeIndex];
        if (existingNode.hasSourceVertexA &&
            incomingNode.hasSourceVertexA &&
            existingNode.sourceVertexIndexA !=
                incomingNode.sourceVertexIndexA) {
            errorMessage = "correspondence merge failure: different chainA vertices collapse to one node";
            return false;
        }
        if (existingNode.hasSourceVertexB &&
            incomingNode.hasSourceVertexB &&
            existingNode.sourceVertexIndexB !=
                incomingNode.sourceVertexIndexB) {
            errorMessage = "correspondence merge failure: different chainB vertices collapse to one node";
            return false;
        }

        const bool keepExistingGeometry =
            existingNode.snappedPair ||
            (!incomingNode.snappedPair &&
                existingNode.projectionDistance <=
                    incomingNode.projectionDistance);
        if (incomingNode.hasSourceVertexA &&
            !existingNode.hasSourceVertexA) {
            existingNode.stationA = incomingNode.stationA;
            existingNode.pointA = incomingNode.pointA;
            existingNode.sourceVertexIndexA =
                incomingNode.sourceVertexIndexA;
            existingNode.hasSourceVertexA = true;
        }
        else if (!existingNode.hasSourceVertexA &&
                 !incomingNode.hasSourceVertexA &&
                 !keepExistingGeometry) {
            existingNode.stationA = incomingNode.stationA;
            existingNode.pointA = incomingNode.pointA;
        }

        if (incomingNode.hasSourceVertexB &&
            !existingNode.hasSourceVertexB) {
            existingNode.stationB = incomingNode.stationB;
            existingNode.pointB = incomingNode.pointB;
            existingNode.sourceVertexIndexB =
                incomingNode.sourceVertexIndexB;
            existingNode.hasSourceVertexB = true;
        }
        else if (!existingNode.hasSourceVertexB &&
                 !incomingNode.hasSourceVertexB &&
                 !keepExistingGeometry) {
            existingNode.stationB = incomingNode.stationB;
            existingNode.pointB = incomingNode.pointB;
        }

        existingNode.projectionDistance = std::min(
            existingNode.projectionDistance,
            incomingNode.projectionDistance);
        existingNode.snappedPair =
            existingNode.snappedPair || incomingNode.snappedPair;
    }

    std::sort(
        mergedNodes.begin(),
        mergedNodes.end(),
        [](const ProjectionCorrespondenceNode& lhs,
           const ProjectionCorrespondenceNode& rhs) {
            if (lhs.stationA != rhs.stationA) {
                return lhs.stationA < rhs.stationA;
            }
            return lhs.stationB < rhs.stationB;
        });

    if (mergedNodes.size() < 2 ||
        !mergedNodes.front().endpoint ||
        !mergedNodes.back().endpoint ||
        std::abs(mergedNodes.front().stationA) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.front().stationB) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.back().stationA - 1.0) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.back().stationB - 1.0) >
            kStationMergeTolerance) {
        errorMessage = "correspondence validation failure: endpoint nodes are invalid";
        return false;
    }
    mergedNodes.front().stationA = 0.0;
    mergedNodes.front().stationB = 0.0;
    mergedNodes.back().stationA = 1.0;
    mergedNodes.back().stationB = 1.0;

    // 对应关系必须在两条 chain 上同时前进；投影回退不能靠重新排序掩盖。
    for (std::size_t nodeIndex = 1;
         nodeIndex < mergedNodes.size();
         ++nodeIndex) {
        const ProjectionCorrespondenceNode& previous =
            mergedNodes[nodeIndex - 1];
        const ProjectionCorrespondenceNode& current =
            mergedNodes[nodeIndex];
        if (current.stationA - previous.stationA <=
                stationMergeToleranceA ||
            current.stationB - previous.stationB <=
                stationMergeToleranceB) {
            errorMessage = "correspondence validation failure: nodes are not strictly monotonic on both chains";
            return false;
        }
    }

    std::vector<int> sourceVertexUseCountA(
        chainA.vertices.size(),
        0);
    std::vector<int> sourceVertexUseCountB(
        chainB.vertices.size(),
        0);
    for (const ProjectionCorrespondenceNode& node : mergedNodes) {
        if (node.pointA.Distance(node.pointB) >
            params_.maxGapWidth + params_.midlineRepairTolerance) {
            errorMessage = "correspondence validation failure: node exceeds allowed gap width";
            return false;
        }
        if (node.hasSourceVertexA) {
            if (node.sourceVertexIndexA == 0 ||
                node.sourceVertexIndexA + 1 >=
                    chainA.vertices.size()) {
                errorMessage = "correspondence validation failure: chainA source vertex index is invalid";
                return false;
            }
            ++sourceVertexUseCountA[node.sourceVertexIndexA];
        }
        if (node.hasSourceVertexB) {
            if (node.sourceVertexIndexB == 0 ||
                node.sourceVertexIndexB + 1 >=
                    chainB.vertices.size()) {
                errorMessage = "correspondence validation failure: chainB source vertex index is invalid";
                return false;
            }
            ++sourceVertexUseCountB[node.sourceVertexIndexB];
        }
    }

    for (std::size_t vertexIndexA = 1;
         vertexIndexA + 1 < sourceVertexUseCountA.size();
         ++vertexIndexA) {
        if (sourceVertexUseCountA[vertexIndexA] != 1) {
            errorMessage = "correspondence validation failure: a chainA internal vertex is not mapped exactly once";
            return false;
        }
    }
    for (std::size_t vertexIndexB = 1;
         vertexIndexB + 1 < sourceVertexUseCountB.size();
         ++vertexIndexB) {
        if (sourceVertexUseCountB[vertexIndexB] != 1) {
            errorMessage = "correspondence validation failure: a chainB internal vertex is not mapped exactly once";
            return false;
        }
    }

    nodes = mergedNodes;
    return true;
}

bool HealClosedNarrowGap::buildProjectionTargetVertices(
    const ClosedGapCandidate& candidate,
    const std::vector<ProjectionCorrespondenceNode>& nodes,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    if (nodes.size() < 2) {
        errorMessage = "target vertex construction failure: correspondence nodes are incomplete";
        return false;
    }

    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    midlineData.chain.vertices.clear();
    midlineData.chain.vertices.reserve(nodes.size());
    midlineData.targetVertexIndicesForChainA.assign(
        candidate.chainA.vertices.size(),
        invalidIndex);
    midlineData.targetVertexIndicesForChainB.assign(
        candidate.chainB.vertices.size(),
        invalidIndex);

    try {
        BRep_Builder builder;
        for (std::size_t nodeIndex = 0;
             nodeIndex < nodes.size();
             ++nodeIndex) {
            if (nodeIndex == 0) {
                midlineData.chain.vertices.push_back(
                    candidate.chainA.vertices.front());
                midlineData.targetVertexIndicesForChainA.front() = 0;
                midlineData.targetVertexIndicesForChainB.front() = 0;
                continue;
            }
            if (nodeIndex + 1 == nodes.size()) {
                midlineData.chain.vertices.push_back(
                    candidate.chainA.vertices.back());
                midlineData.targetVertexIndicesForChainA.back() =
                    nodeIndex;
                midlineData.targetVertexIndicesForChainB.back() =
                    nodeIndex;
                continue;
            }

            const ProjectionCorrespondenceNode& node =
                nodes[nodeIndex];
            // snapped pair 使用两个真实顶点；普通节点使用源顶点和对侧投影点，
            // 两种情况都取二者中点，避免把整条缝隙强行压到任意一侧。
            const gp_Pnt targetPoint(
                0.5 * (node.pointA.X() + node.pointB.X()),
                0.5 * (node.pointA.Y() + node.pointB.Y()),
                0.5 * (node.pointA.Z() + node.pointB.Z()));
            if (!isFinitePoint(targetPoint)) {
                errorMessage = "target vertex construction failure: midpoint is not finite";
                return false;
            }

            double vertexTolerance = std::max({
                Precision::Confusion(),
                params_.midlineRepairTolerance,
                node.pointA.Distance(targetPoint)
                    + Precision::Confusion(),
                node.pointB.Distance(targetPoint)
                    + Precision::Confusion()
            });
            if (node.hasSourceVertexA) {
                const TopoDS_Vertex& sourceVertex =
                    candidate.chainA.vertices[
                        node.sourceVertexIndexA];
                vertexTolerance = std::max({
                    vertexTolerance,
                    BRep_Tool::Tolerance(sourceVertex),
                    BRep_Tool::Pnt(sourceVertex).Distance(targetPoint)
                        + Precision::Confusion()
                });
            }
            if (node.hasSourceVertexB) {
                const TopoDS_Vertex& sourceVertex =
                    candidate.chainB.vertices[
                        node.sourceVertexIndexB];
                vertexTolerance = std::max({
                    vertexTolerance,
                    BRep_Tool::Tolerance(sourceVertex),
                    BRep_Tool::Pnt(sourceVertex).Distance(targetPoint)
                        + Precision::Confusion()
                });
            }

            TopoDS_Vertex targetVertex;
            builder.MakeVertex(
                targetVertex,
                targetPoint,
                vertexTolerance);
            if (targetVertex.IsNull()) {
                errorMessage = "target vertex construction failure: failed to create midpoint vertex";
                return false;
            }
            midlineData.chain.vertices.push_back(targetVertex);

            if (node.hasSourceVertexA) {
                std::size_t& mappedIndex =
                    midlineData.targetVertexIndicesForChainA[
                        node.sourceVertexIndexA];
                if (mappedIndex != invalidIndex) {
                    errorMessage = "target vertex construction failure: chainA vertex maps more than once";
                    return false;
                }
                mappedIndex = nodeIndex;
            }
            if (node.hasSourceVertexB) {
                std::size_t& mappedIndex =
                    midlineData.targetVertexIndicesForChainB[
                        node.sourceVertexIndexB];
                if (mappedIndex != invalidIndex) {
                    errorMessage = "target vertex construction failure: chainB vertex maps more than once";
                    return false;
                }
                mappedIndex = nodeIndex;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "target vertex construction failure: "
            + failureText(failure);
        return false;
    }

    for (std::size_t vertexIndexA = 0;
         vertexIndexA <
            midlineData.targetVertexIndicesForChainA.size();
         ++vertexIndexA) {
        if (midlineData.targetVertexIndicesForChainA[vertexIndexA] ==
            invalidIndex) {
            errorMessage = "target vertex construction failure: chainA vertex mapping is incomplete";
            return false;
        }
    }
    for (std::size_t vertexIndexB = 0;
         vertexIndexB <
            midlineData.targetVertexIndicesForChainB.size();
         ++vertexIndexB) {
        if (midlineData.targetVertexIndicesForChainB[vertexIndexB] ==
            invalidIndex) {
            errorMessage = "target vertex construction failure: chainB vertex mapping is incomplete";
            return false;
        }
    }
    if (midlineData.chain.vertices.size() != nodes.size()) {
        errorMessage = "target vertex construction failure: target vertex count is inconsistent";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::buildProjectionTargetEdges(
    const ClosedGapCandidate& candidate,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    const std::vector<ProjectionCorrespondenceNode>& nodes,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    if (nodes.size() < 2 ||
        midlineData.chain.vertices.size() != nodes.size()) {
        errorMessage = "target edge construction failure: target vertices are inconsistent";
        return false;
    }

    midlineData.chain.edges.clear();
    midlineData.chain.length = 0.0;
    midlineData.chain.edges.reserve(nodes.size() - 1);
    try {
        for (std::size_t targetEdgeIndex = 0;
             targetEdgeIndex + 1 < nodes.size();
             ++targetEdgeIndex) {
            const ProjectionCorrespondenceNode& startNode =
                nodes[targetEdgeIndex];
            const ProjectionCorrespondenceNode& endNode =
                nodes[targetEdgeIndex + 1];
            const double stationSpanA =
                endNode.stationA - startNode.stationA;
            const double stationSpanB =
                endNode.stationB - startNode.stationB;
            if (!std::isfinite(stationSpanA) ||
                !std::isfinite(stationSpanB) ||
                stationSpanA <= 0.0 ||
                stationSpanB <= 0.0) {
                errorMessage = "target edge construction failure: correspondence interval is invalid";
                return false;
            }

            const double segmentLengthA =
                chainAStations.totalLength * stationSpanA;
            const double segmentLengthB =
                chainBStations.totalLength * stationSpanB;
            const int sampleCount = computeMidlineSampleCount(
                segmentLengthA,
                segmentLengthB);
            if (sampleCount < 2) {
                errorMessage = "target edge construction failure: sample count is invalid";
                return false;
            }

            const gp_Pnt startPoint = BRep_Tool::Pnt(
                midlineData.chain.vertices[targetEdgeIndex]);
            const gp_Pnt endPoint = BRep_Tool::Pnt(
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            if (!isFinitePoint(startPoint) ||
                !isFinitePoint(endPoint) ||
                startPoint.Distance(endPoint) <=
                    Precision::Confusion()) {
                errorMessage = "target edge construction failure: target endpoints are invalid";
                return false;
            }

            std::vector<gp_Pnt> guidancePoints;
            guidancePoints.reserve(
                static_cast<std::size_t>(sampleCount));
            guidancePoints.push_back(startPoint);
            for (int sampleIndex = 1;
                 sampleIndex + 1 < sampleCount;
                 ++sampleIndex) {
                const double local =
                    static_cast<double>(sampleIndex)
                    / static_cast<double>(sampleCount - 1);
                const double stationA =
                    startNode.stationA + stationSpanA * local;
                const double stationB =
                    startNode.stationB + stationSpanB * local;

                gp_Pnt pointA;
                gp_Pnt pointB;
                if (!sampleChainPoint(
                        candidate.chainA,
                        stationA,
                        pointA) ||
                    !sampleChainPoint(
                        candidate.chainB,
                        stationB,
                        pointB)) {
                    errorMessage = "target edge construction failure: failed to sample correspondence interval";
                    return false;
                }
                if (!isFinitePoint(pointA) ||
                    !isFinitePoint(pointB) ||
                    pointA.Distance(pointB) >
                        params_.maxGapWidth
                            + params_.midlineRepairTolerance) {
                    errorMessage = "target edge construction failure: interval samples exceed allowed gap width";
                    return false;
                }

                // A、B 使用各自站点插值，表达投影对应而不是错误的相同 t 对应。
                const gp_Pnt midpoint(
                    0.5 * (pointA.X() + pointB.X()),
                    0.5 * (pointA.Y() + pointB.Y()),
                    0.5 * (pointA.Z() + pointB.Z()));
                if (!isFinitePoint(midpoint)) {
                    errorMessage = "target edge construction failure: sample midpoint is invalid";
                    return false;
                }
                if (guidancePoints.back().Distance(midpoint) >
                        Precision::Confusion() &&
                    midpoint.Distance(endPoint) >
                        Precision::Confusion()) {
                    guidancePoints.push_back(midpoint);
                }
            }

            if (guidancePoints.back().Distance(endPoint) <=
                Precision::Confusion()) {
                errorMessage = "target edge construction failure: guidance points collapse at endpoint";
                return false;
            }
            guidancePoints.push_back(endPoint);

            Handle(TColgp_HArray1OfPnt) interpolationPoints =
                new TColgp_HArray1OfPnt(
                    1,
                    static_cast<Standard_Integer>(
                        guidancePoints.size()));
            for (std::size_t pointIndex = 0;
                 pointIndex < guidancePoints.size();
                 ++pointIndex) {
                interpolationPoints->SetValue(
                    static_cast<Standard_Integer>(pointIndex + 1),
                    guidancePoints[pointIndex]);
            }

            GeomAPI_Interpolate interpolator(
                interpolationPoints,
                Standard_False,
                Precision::Confusion());
            interpolator.Perform();
            if (!interpolator.IsDone() ||
                interpolator.Curve().IsNull()) {
                errorMessage = "target edge construction failure: interpolation failed";
                return false;
            }

            BRepBuilderAPI_MakeEdge edgeMaker(
                interpolator.Curve(),
                midlineData.chain.vertices[targetEdgeIndex],
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "target edge construction failure: edge builder failed";
                return false;
            }

            TopoDS_Edge targetEdge = orientEdgeFromTo(
                edgeMaker.Edge(),
                midlineData.chain.vertices[targetEdgeIndex],
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            TopoDS_Vertex targetStart;
            TopoDS_Vertex targetEnd;
            if (targetEdge.IsNull() ||
                !getEdgeVertices(
                    targetEdge,
                    targetStart,
                    targetEnd) ||
                !targetStart.IsSame(
                    midlineData.chain.vertices[targetEdgeIndex]) ||
                !targetEnd.IsSame(
                    midlineData.chain.vertices[
                        targetEdgeIndex + 1])) {
                errorMessage = "target edge construction failure: canonical edge orientation is invalid";
                return false;
            }

            const double targetEdgeLength =
                GeometryUtils::ComputeEdgeLength(targetEdge);
            if (!std::isfinite(targetEdgeLength) ||
                targetEdgeLength <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: canonical edge is degenerate";
                return false;
            }

            // canonical target edge 只创建一次，后续由 A、B 两侧 face 共同引用。
            midlineData.chain.edges.push_back(targetEdge);
            midlineData.chain.length += targetEdgeLength;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "target edge construction failure: "
            + failureText(failure);
        return false;
    }

    if (midlineData.chain.edges.size() + 1 !=
            midlineData.chain.vertices.size() ||
        !std::isfinite(midlineData.chain.length) ||
        midlineData.chain.length <= Precision::Confusion()) {
        errorMessage = "target edge construction failure: canonical target chain is invalid";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::mergeProjectionNodesForShortTargetEdges(
    std::vector<ProjectionCorrespondenceNode>& nodes,
    const MidlineChainData& midlineData,
    bool& nodesMerged,
    std::string& errorMessage) const
{
    nodesMerged = false;
    if (nodes.size() < 2 ||
        midlineData.chain.vertices.size() != nodes.size() ||
        midlineData.chain.edges.size() + 1 != nodes.size()) {
        errorMessage = "short target edge merge failure: target topology is inconsistent";
        return false;
    }

    try {
        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.chain.edges.size();
             ++edgeIndex) {
            const TopoDS_Edge& targetEdge =
                midlineData.chain.edges[edgeIndex];
            const TopoDS_Vertex& startVertex =
                midlineData.chain.vertices[edgeIndex];
            const TopoDS_Vertex& endVertex =
                midlineData.chain.vertices[edgeIndex + 1];
            const double edgeLength =
                GeometryUtils::ComputeEdgeLength(targetEdge);
            const double endpointTolerance =
                std::max(Precision::Confusion(),
                    BRep_Tool::Tolerance(startVertex)
                        + BRep_Tool::Tolerance(endVertex));
            const bool mayDegenerate =
                BRep_Tool::Degenerated(targetEdge) == Standard_True ||
                !std::isfinite(edgeLength) ||
                edgeLength <= endpointTolerance;
            if (!mayDegenerate) {
                continue;
            }

            // 端点节点代表两侧 chain 的公共尖点，不能为了消除短边而吞并内部拓扑。
            if (edgeIndex == 0 ||
                edgeIndex + 2 == nodes.size()) {
                std::ostringstream os;
                os << "short target edge merge failure: near-degenerate edge "
                   << edgeIndex
                   << " touches a chain endpoint; length=" << edgeLength
                   << ", endpointTolerance=" << endpointTolerance;
                errorMessage = os.str();
                return false;
            }

            const ProjectionCorrespondenceNode& firstNode =
                nodes[edgeIndex];
            const ProjectionCorrespondenceNode& secondNode =
                nodes[edgeIndex + 1];
            if ((firstNode.hasSourceVertexA &&
                 secondNode.hasSourceVertexA &&
                 firstNode.sourceVertexIndexA !=
                     secondNode.sourceVertexIndexA) ||
                (firstNode.hasSourceVertexB &&
                 secondNode.hasSourceVertexB &&
                 firstNode.sourceVertexIndexB !=
                     secondNode.sourceVertexIndexB)) {
                std::ostringstream os;
                os << "short target edge merge failure: edge "
                   << edgeIndex
                   << " cannot be removed without collapsing distinct "
                      "source vertices on the same chain; length="
                   << edgeLength
                   << ", endpointTolerance=" << endpointTolerance;
                errorMessage = os.str();
                return false;
            }

            const bool preferFirst =
                firstNode.snappedPair ||
                (!secondNode.snappedPair &&
                 firstNode.projectionDistance <=
                     secondNode.projectionDistance);
            ProjectionCorrespondenceNode mergedNode =
                preferFirst ? firstNode : secondNode;

            // 每一侧都优先保留该侧真实 source vertex 的站点和坐标。
            // 对互补的 A-only/B-only 节点，这会恢复成一个完整的 A/B 对应节点。
            if (firstNode.hasSourceVertexA) {
                mergedNode.stationA = firstNode.stationA;
                mergedNode.pointA = firstNode.pointA;
                mergedNode.sourceVertexIndexA =
                    firstNode.sourceVertexIndexA;
                mergedNode.hasSourceVertexA = true;
            }
            else if (secondNode.hasSourceVertexA) {
                mergedNode.stationA = secondNode.stationA;
                mergedNode.pointA = secondNode.pointA;
                mergedNode.sourceVertexIndexA =
                    secondNode.sourceVertexIndexA;
                mergedNode.hasSourceVertexA = true;
            }

            if (firstNode.hasSourceVertexB) {
                mergedNode.stationB = firstNode.stationB;
                mergedNode.pointB = firstNode.pointB;
                mergedNode.sourceVertexIndexB =
                    firstNode.sourceVertexIndexB;
                mergedNode.hasSourceVertexB = true;
            }
            else if (secondNode.hasSourceVertexB) {
                mergedNode.stationB = secondNode.stationB;
                mergedNode.pointB = secondNode.pointB;
                mergedNode.sourceVertexIndexB =
                    secondNode.sourceVertexIndexB;
                mergedNode.hasSourceVertexB = true;
            }

            mergedNode.endpoint = false;
            mergedNode.snappedPair =
                firstNode.snappedPair || secondNode.snappedPair;
            mergedNode.mergedForShortEdge = true;
            mergedNode.projectionDistance =
                mergedNode.pointA.Distance(mergedNode.pointB);
            if (!std::isfinite(mergedNode.projectionDistance) ||
                mergedNode.projectionDistance >
                    params_.maxGapWidth
                        + params_.midlineRepairTolerance) {
                errorMessage = "short target edge merge failure: merged node exceeds allowed gap width";
                return false;
            }

            const ProjectionCorrespondenceNode& previousNode =
                nodes[edgeIndex - 1];
            const ProjectionCorrespondenceNode& nextNode =
                nodes[edgeIndex + 2];
            if (mergedNode.stationA - previousNode.stationA <=
                    kStationMergeTolerance ||
                mergedNode.stationB - previousNode.stationB <=
                    kStationMergeTolerance ||
                nextNode.stationA - mergedNode.stationA <=
                    kStationMergeTolerance ||
                nextNode.stationB - mergedNode.stationB <=
                    kStationMergeTolerance) {
                errorMessage = "short target edge merge failure: merged node breaks correspondence monotonicity";
                return false;
            }

            nodes[edgeIndex] = mergedNode;
            nodes.erase(nodes.begin()
                + static_cast<std::ptrdiff_t>(edgeIndex + 1));
            nodesMerged = true;
            return true;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "short target edge merge failure: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::buildSourceToTargetEdgeMapping(
    const std::vector<double>& sourceVertexStations,
    const std::vector<double>& targetStations,
    const char* chainName,
    std::vector<std::vector<std::size_t>>& mapping,
    std::string& errorMessage) const
{
    mapping.clear();
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (sourceVertexStations.size() < 2 || targetStations.size() < 2) {
        errorMessage = "mapping failure: " + name
            + " station sequence is incomplete";
        return false;
    }

    mapping.resize(sourceVertexStations.size() - 1);
    for (std::size_t targetEdgeIndex = 0;
         targetEdgeIndex + 1 < targetStations.size();
         ++targetEdgeIndex) {
        const double intervalStart = targetStations[targetEdgeIndex];
        const double intervalEnd = targetStations[targetEdgeIndex + 1];
        const double midpoint = 0.5 * (intervalStart + intervalEnd);
        if (!std::isfinite(midpoint) ||
            midpoint <= intervalStart ||
            midpoint >= intervalEnd) {
            errorMessage = "mapping failure: target edge interval is invalid";
            return false;
        }

        const auto upper = std::upper_bound(
            sourceVertexStations.begin(),
            sourceVertexStations.end(),
            midpoint);
        if (upper == sourceVertexStations.begin() ||
            upper == sourceVertexStations.end()) {
            errorMessage = "mapping failure: " + name
                + " cannot locate the owner source edge for a target edge";
            return false;
        }

        const std::size_t sourceEdgeIndex = static_cast<std::size_t>(
            std::distance(sourceVertexStations.begin(), upper) - 1);
        if (sourceEdgeIndex >= mapping.size() ||
            midpoint <= sourceVertexStations[sourceEdgeIndex] ||
            midpoint >= sourceVertexStations[sourceEdgeIndex + 1]) {
            errorMessage = "mapping failure: " + name
                + " target edge midpoint is outside its source interval";
            return false;
        }
        mapping[sourceEdgeIndex].push_back(targetEdgeIndex);
    }

    return true;
}

bool HealClosedNarrowGap::validateReplacementMapping(
    const std::vector<std::vector<std::size_t>>& mapping,
    std::size_t sourceEdgeCount,
    std::size_t targetEdgeCount,
    const char* chainName,
    std::string& errorMessage) const
{
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (mapping.size() != sourceEdgeCount) {
        errorMessage = "mapping failure: " + name
            + " replacement mapping count does not match source edges";
        return false;
    }
    if (targetEdgeCount == 0) {
        errorMessage = "mapping failure: canonical target edge chain is empty";
        return false;
    }

    std::vector<int> targetOwnership(targetEdgeCount, 0);
    std::size_t nextExpectedTargetIndex = 0;
    for (std::size_t sourceEdgeIndex = 0;
         sourceEdgeIndex < mapping.size();
         ++sourceEdgeIndex) {
        const std::vector<std::size_t>& replacements =
            mapping[sourceEdgeIndex];
        if (replacements.empty()) {
            errorMessage = "mapping failure: " + name
                + " source edge has no target edge";
            return false;
        }

        for (std::size_t replacementIndex = 0;
             replacementIndex < replacements.size();
             ++replacementIndex) {
            const std::size_t targetEdgeIndex =
                replacements[replacementIndex];
            if (targetEdgeIndex >= targetEdgeCount) {
                errorMessage = "mapping failure: " + name
                    + " target edge index is out of range";
                return false;
            }
            if (targetEdgeIndex != nextExpectedTargetIndex) {
                errorMessage = "mapping failure: " + name
                    + " target edges are not ordered by source edge";
                return false;
            }
            if (replacementIndex > 0 &&
                targetEdgeIndex != replacements[replacementIndex - 1] + 1) {
                errorMessage = "mapping failure: " + name
                    + " replacement indices are not continuous and increasing";
                return false;
            }
            ++targetOwnership[targetEdgeIndex];
            ++nextExpectedTargetIndex;
        }
    }

    if (nextExpectedTargetIndex != targetEdgeCount) {
        errorMessage = "mapping failure: " + name
            + " does not cover the full canonical target chain";
        return false;
    }
    for (std::size_t targetEdgeIndex = 0;
         targetEdgeIndex < targetOwnership.size();
         ++targetEdgeIndex) {
        if (targetOwnership[targetEdgeIndex] != 1) {
            errorMessage = "mapping failure: " + name
                + " does not cover every target edge exactly once";
            return false;
        }
    }
    return true;
}

bool HealClosedNarrowGap::buildMidlineChainByProjectionCorrespondence(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    midlineData.constructionMode =
        MidlineConstructionMode::ProjectionCorrespondence;

    ChainStationData chainAStations;
    ChainStationData chainBStations;
    if (!buildNormalizedVertexStations(
            candidate.chainA,
            "chainA",
            chainAStations,
            errorMessage) ||
        !buildNormalizedVertexStations(
            candidate.chainB,
            "chainB",
            chainBStations,
            errorMessage)) {
        return false;
    }

    std::vector<ProjectionCorrespondenceNode> nodes;
    if (!buildProjectionCorrespondenceNodes(
            candidate,
            chainAStations,
            chainBStations,
            nodes,
            errorMessage)) {
        return false;
    }

    const std::size_t initialNodeCount = nodes.size();
    while (true) {
        if (!buildProjectionTargetVertices(
                candidate,
                nodes,
                midlineData,
                errorMessage) ||
            !buildProjectionTargetEdges(
                candidate,
                chainAStations,
                chainBStations,
                nodes,
                midlineData,
                errorMessage)) {
            return false;
        }

        bool nodesMerged = false;
        if (!mergeProjectionNodesForShortTargetEdges(
                nodes,
                midlineData,
                nodesMerged,
                errorMessage)) {
            return false;
        }
        if (!nodesMerged) {
            break;
        }

        ++midlineData.shortEdgeMergeCount;
        if (nodes.size() < 2 ||
            midlineData.shortEdgeMergeCount >
                static_cast<int>(initialNodeCount - 2)) {
            errorMessage = "short target edge merge failure: merge iteration limit exceeded";
            return false;
        }
    }
    midlineData.correspondenceNodes = nodes;

    std::vector<double> targetStationsA;
    std::vector<double> targetStationsB;
    targetStationsA.reserve(nodes.size());
    targetStationsB.reserve(nodes.size());
    for (const ProjectionCorrespondenceNode& node : nodes) {
        targetStationsA.push_back(node.stationA);
        targetStationsB.push_back(node.stationB);
    }

    // 一条 source edge 可以跨越多个对应区间，但每个 target edge
    // 在同一侧只能归属于一条 source edge。
    if (!buildSourceToTargetEdgeMapping(
            chainAStations.vertexStations,
            targetStationsA,
            "chainA",
            midlineData.replacementsForChainA,
            errorMessage) ||
        !buildSourceToTargetEdgeMapping(
            chainBStations.vertexStations,
            targetStationsB,
            "chainB",
            midlineData.replacementsForChainB,
            errorMessage)) {
        return false;
    }

    if (!validateReplacementMapping(
            midlineData.replacementsForChainA,
            candidate.chainA.edges.size(),
            midlineData.chain.edges.size(),
            "chainA",
            errorMessage) ||
        !validateReplacementMapping(
            midlineData.replacementsForChainB,
            candidate.chainB.edges.size(),
            midlineData.chain.edges.size(),
            "chainB",
            errorMessage)) {
        return false;
    }

    return true;
}

const char* HealClosedNarrowGap::midlineConstructionModeName(
    MidlineConstructionMode mode)
{
    switch (mode) {
    case MidlineConstructionMode::Auto:
        return "Auto";
    case MidlineConstructionMode::LegacyOneToOne:
        return "LegacyOneToOne";
    case MidlineConstructionMode::ProjectionCorrespondence:
        return "ProjectionCorrespondence";
    default:
        return "Unknown";
    }
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
        const bool asymmetricSourceChains =
            candidate.chainA.edges.size() != candidate.chainB.edges.size();
        if (asymmetricSourceChains) {
            sideError =
                "candidate uses asymmetric source chains; "
                "midline collapse is disabled or failed; "
                "legacy side replacement only supports one-to-one source edge replacement";
        }
        else if (tryRepairBySideReplacement(
                     currentShape,
                     candidate,
                     repairedShape,
                     sideError)) {
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

    MidlineConstructionMode resolvedMode =
        MidlineConstructionMode::Auto;
    std::string modeError;
    if (!resolveMidlineConstructionMode(
            candidate,
            resolvedMode,
            modeError)) {
        appendMessage(
            "Midline construction mode: "
            + std::string(midlineConstructionModeName(
                params_.midlineConstructionMode)));
        errorMessage = modeError;
        return false;
    }
    appendMessage(
        "Midline construction mode: "
        + std::string(midlineConstructionModeName(resolvedMode)));

    MidlineChainData midlineData;
    if (!buildMidlineChain(candidate, midlineData, errorMessage)) {
        return false;
    }

    if (midlineData.constructionMode ==
        MidlineConstructionMode::ProjectionCorrespondence) {
        std::ostringstream summary;
        summary << "Projection correspondence midline topology: chainAEdges="
                << candidate.chainA.edges.size()
                << ", chainBEdges=" << candidate.chainB.edges.size()
                << ", chainAVertices=" << candidate.chainA.vertices.size()
                << ", chainBVertices=" << candidate.chainB.vertices.size()
                << ", correspondenceNodeCount="
                << midlineData.correspondenceNodes.size()
                << ", shortEdgeNodeMerges="
                << midlineData.shortEdgeMergeCount
                << ", targetVertices=" << midlineData.chain.vertices.size()
                << ", targetEdges=" << midlineData.chain.edges.size();
        appendMessage(summary.str());

        for (std::size_t nodeIndex = 0;
             nodeIndex < midlineData.correspondenceNodes.size();
             ++nodeIndex) {
            const ProjectionCorrespondenceNode& node =
                midlineData.correspondenceNodes[nodeIndex];
            std::ostringstream nodeLog;
            nodeLog << "correspondence node " << nodeIndex
                    << ": stationA=" << node.stationA
                    << ", stationB=" << node.stationB
                    << ", endpoint="
                    << (node.endpoint ? "true" : "false")
                    << ", sourceVertexA=";
            if (node.hasSourceVertexA) {
                nodeLog << node.sourceVertexIndexA;
            }
            else {
                nodeLog << "none";
            }
            nodeLog << ", sourceVertexB=";
            if (node.hasSourceVertexB) {
                nodeLog << node.sourceVertexIndexB;
            }
            else {
                nodeLog << "none";
            }
            nodeLog << ", snappedPair="
                    << (node.snappedPair ? "true" : "false")
                    << ", mergedForShortEdge="
                    << (node.mergedForShortEdge ? "true" : "false")
                    << ", projectionDistance="
                    << node.projectionDistance;
            appendMessage(nodeLog.str());
        }

        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.replacementsForChainA.size();
             ++edgeIndex) {
            std::ostringstream mapping;
            mapping << "chainA edge " << edgeIndex
                    << " -> target edges "
                    << formatIndexList(
                        midlineData.replacementsForChainA[edgeIndex]);
            appendMessage(mapping.str());
        }
        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.replacementsForChainB.size();
             ++edgeIndex) {
            std::ostringstream mapping;
            mapping << "chainB edge " << edgeIndex
                    << " -> target edges "
                    << formatIndexList(
                        midlineData.replacementsForChainB[edgeIndex]);
            appendMessage(mapping.str());
        }
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
        errorMessage =
            "candidate uses asymmetric source chains; "
            "midline collapse is disabled or failed; "
            "legacy side replacement only supports one-to-one source edge replacement";
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
    if (midlineData.chain.vertices.size() !=
        midlineData.chain.edges.size() + 1) {
        errorMessage = "midline target vertex and edge counts are inconsistent";
        return false;
    }
    if (!validateReplacementMapping(
            midlineData.replacementsForChainA,
            candidate.chainA.edges.size(),
            midlineData.chain.edges.size(),
            "chainA",
            errorMessage) ||
        !validateReplacementMapping(
            midlineData.replacementsForChainB,
            candidate.chainB.edges.size(),
            midlineData.chain.edges.size(),
            "chainB",
            errorMessage)) {
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
            for (std::size_t replacementIndex = 0;
                 replacementIndex < replacements[edgeIndex].size();
                 ++replacementIndex) {
                const std::size_t targetEdgeIndex =
                    replacements[edgeIndex][replacementIndex];
                if (targetEdgeIndex >= midlineData.chain.edges.size()) {
                    errorMessage = std::string(chainName)
                        + " replacement target edge index is out of range";
                    return false;
                }
                if (replacementIndex > 0 &&
                    targetEdgeIndex !=
                        replacements[edgeIndex][replacementIndex - 1] + 1) {
                    errorMessage = std::string(chainName)
                        + " replacement target edge indices are not continuous";
                    return false;
                }
            }

            FaceReplacementPlan* plan = findOrAddFacePlan(ownerFace);
            for (const EdgeReplacementTask& existingTask : plan->edgeTasks) {
                if (existingTask.sourceEdge.IsSame(sourceChain.edges[edgeIndex])) {
                    errorMessage = "the same source edge appears more than once in midline replacement tasks";
                    return false;
                }
                for (const std::size_t existingTargetIndex :
                     existingTask.replacementEdgeIndices) {
                    if (std::find(
                            replacements[edgeIndex].begin(),
                            replacements[edgeIndex].end(),
                            existingTargetIndex)
                        != replacements[edgeIndex].end()) {
                        errorMessage = "the same canonical target edge would be inserted "
                            "more than once into one face";
                        return false;
                    }
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

        const std::size_t invalidTargetVertexIndex =
            std::numeric_limits<std::size_t>::max();
        auto registerChainVertexReplacements =
            [&](const EdgeChain& sourceChain,
                const std::vector<std::size_t>& targetVertexIndices,
                const char* chainName) {
                if (targetVertexIndices.size() !=
                    sourceChain.vertices.size()) {
                    errorMessage = std::string(
                        "vertex replacement failure: ")
                        + chainName
                        + " source vertex mapping count is inconsistent";
                    return false;
                }
                if (targetVertexIndices.empty() ||
                    targetVertexIndices.front() != 0 ||
                    targetVertexIndices.back() ==
                        invalidTargetVertexIndex ||
                    targetVertexIndices.back() !=
                        midlineData.chain.vertices.size() - 1) {
                    errorMessage = std::string(
                        "vertex replacement failure: ")
                        + chainName
                        + " endpoint mapping is inconsistent";
                    return false;
                }

                // A、B 两侧内部顶点数量可能不同，分别按各自映射替换，
                // 不能再假设两侧 source vertex 与 target vertex 使用相同下标。
                for (std::size_t sourceVertexIndex = 1;
                     sourceVertexIndex + 1 < sourceChain.vertices.size();
                     ++sourceVertexIndex) {
                    const std::size_t targetVertexIndex =
                        targetVertexIndices[sourceVertexIndex];
                    if (targetVertexIndex ==
                            invalidTargetVertexIndex ||
                        targetVertexIndex == 0 ||
                        targetVertexIndex + 1 >=
                            midlineData.chain.vertices.size()) {
                        errorMessage = std::string(
                            "vertex replacement failure: ")
                            + chainName
                            + " target vertex index is invalid";
                        return false;
                    }

                    if (!registerVertexReplacement(
                            sourceChain.vertices[sourceVertexIndex],
                            midlineData.chain.vertices[
                                targetVertexIndex])) {
                        if (errorMessage.rfind(
                                "vertex replacement failure:",
                                0) != 0) {
                            errorMessage =
                                "vertex replacement failure: "
                                + errorMessage;
                        }
                        return false;
                    }
                }
                return true;
            };

        if (!registerChainVertexReplacements(
                candidate.chainA,
                midlineData.targetVertexIndicesForChainA,
                "chainA") ||
            !registerChainVertexReplacements(
                candidate.chainB,
                midlineData.targetVertexIndicesForChainB,
                "chainB")) {
            return false;
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
                        // source edge 在 wire 中反向使用时，规范目标链的遍历方向也整体反转，
                        // 因此既要倒序插入 replacement edges，也要反转每个 edge use。
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

int HealClosedNarrowGap::computeWidthSampleCount(
    double lengthA,
    double lengthB) const
{
    const int fallbackCount = std::max(2, params_.widthSampleCount);
    if (!params_.useAdaptiveWidthSampleCount ||
        !std::isfinite(params_.targetWidthSampleSpacing) ||
        params_.targetWidthSampleSpacing <= Precision::Confusion()) {
        return fallbackCount;
    }

    if (!std::isfinite(lengthA) ||
        !std::isfinite(lengthB) ||
        lengthA <= Precision::Confusion() ||
        lengthB <= Precision::Confusion()) {
        return fallbackCount;
    }

    const double pairLength = std::max(lengthA, lengthB);
    const double countByLength =
        std::ceil(pairLength / params_.targetWidthSampleSpacing) + 1.0;
    if (!std::isfinite(countByLength) ||
        countByLength >= static_cast<double>(params_.maxWidthSampleCount)) {
        return params_.maxWidthSampleCount;
    }
    if (countByLength <= static_cast<double>(params_.minWidthSampleCount)) {
        return params_.minWidthSampleCount;
    }

    return static_cast<int>(countByLength);
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
