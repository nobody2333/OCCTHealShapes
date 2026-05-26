#include <algorithm>
#include <cmath>
#include <limits>

#include "MergeShortEdge.h"
#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepTools_ReShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <gp_Pnt.hxx>
#include <TopTools_MapOfShape.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Line.hxx>
#include <GeomLib.hxx>
#include <Geom_Circle.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Parabola.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BoundedCurve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Standard_Failure.hxx>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MergeShortEdge::MergeShortEdge(double bboxRatio, double angleThresholdDeg)
    : lengthThreshold_(0.0)
    , angleThresholdDeg_(angleThresholdDeg)
    , bboxRatio_(bboxRatio)
    , localLengthRatio_(0.25)
    , thresholdMode_(bboxRatio > 0.0 ? ShortEdgeThresholdMode::HybridGlobalAndLocal : ShortEdgeThresholdMode::ManualAbsolute)
{
}

void MergeShortEdge::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;

	allEdges_.clear();
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    edgeFaceMap_.Clear();
    resultShape_.Nullify();
    clearVirtualPolyCurves();
    refreshLengthThreshold();
}

void MergeShortEdge::clearVirtualPolyCurves()
{
    virtualPolyCurves_.clear();
    edgePolyCurveMap_.Clear();
}


PolyCurveEdgeUse MergeShortEdge::makePolyCurveEdgeUse(
    const TopoDS_Edge& edge,
    const TopoDS_Vertex& startVertex,
    const TopoDS_Vertex& endVertex) const
{
    PolyCurveEdgeUse use;
    use.edge = edge;
    use.reversed = false;

    if (edge.IsNull() || startVertex.IsNull() || endVertex.IsNull()) {
        return use;
    }

    const TopoDS_Vertex edgeFirst = TopExp::FirstVertex(edge);
    const TopoDS_Vertex edgeLast = TopExp::LastVertex(edge);

    if (edgeFirst.IsSame(startVertex) && edgeLast.IsSame(endVertex)) {
        use.reversed = false;
    }
    else if (edgeFirst.IsSame(endVertex) && edgeLast.IsSame(startVertex)) {
        use.reversed = true;
    }

    return use;
}

VirtualPolyCurve MergeShortEdge::createVirtualPolyCurveFromCandidate(
    const MergeCandidate& candidate,
    int id) const
{
    VirtualPolyCurve curve;
    curve.id = id;

    const TopoDS_Vertex other1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    const TopoDS_Vertex other2 = getOtherVertex(candidate.edge2, candidate.commonVertex);
    if (other1.IsNull() || other2.IsNull()) {
        return curve;
    }

    // 与 createMergedEdgeByCompositeBSpline() 保持同一条 PolyCurve 的方向规则，
    // 但这里只记录 reference edges，不重建 OCCT 边、不删除短边。
    if (candidate.commonVertex.IsSame(TopExp::FirstVertex(candidate.edge2))) {
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge1, other1, candidate.commonVertex));
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge2, candidate.commonVertex, other2));
    }
    else if (candidate.commonVertex.IsSame(TopExp::LastVertex(candidate.edge2))) {
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge2, other2, candidate.commonVertex));
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge1, candidate.commonVertex, other1));
    }
    else {
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge1, other1, candidate.commonVertex));
        curve.edges.push_back(makePolyCurveEdgeUse(candidate.edge2, candidate.commonVertex, other2));
    }

    return curve;
}

int MergeShortEdge::buildVirtualPolyCurvesFromMergeCandidates()
{
    clearVirtualPolyCurves();
    if(mergeCandidates_.empty()) {
        findMergeCandidates();
	}
    int nextId = 1;

    for (const auto& candidate : mergeCandidates_) {
        if (candidate.edge1.IsNull() || candidate.edge2.IsNull()) {
            continue;
        }
        if (edgePolyCurveMap_.IsBound(candidate.edge1) || edgePolyCurveMap_.IsBound(candidate.edge2)) {
            continue;
        }

        VirtualPolyCurve polyCurves = createVirtualPolyCurveFromCandidate(candidate, nextId);
        if (polyCurves.edges.size() < 2) {
            continue;
        }

        for (const auto& edgeUse : polyCurves.edges) {
            if (!edgeUse.edge.IsNull() && !edgePolyCurveMap_.IsBound(edgeUse.edge)) {
                edgePolyCurveMap_.Bind(edgeUse.edge, polyCurves.id);
            }
        }

        virtualPolyCurves_.push_back(polyCurves);
        ++nextId;
    }

    return static_cast<int>(virtualPolyCurves_.size());
}

void MergeShortEdge::buildTopologyMap()
{
    vertexEdgeMap_.Clear();
    edgeFaceMap_.Clear();
    if (inputShape_.IsNull()) {
        return;
    }

    TopologyUtils::BuildVertexEdgeMap(inputShape_, vertexEdgeMap_);
    TopExp::MapShapesAndAncestors(inputShape_, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap_);
}

bool MergeShortEdge::refreshLengthThreshold()
{
    // 统一入口：后续 setInput / performMerge 只调用这里，确保阈值行为可预期。
    switch (thresholdMode_) {
    case ShortEdgeThresholdMode::ManualAbsolute:
        return lengthThreshold_ > 0.0;
    case ShortEdgeThresholdMode::BoundingBoxRatio:
        return setLengthThresholdByBoundingBox(bboxRatio_);
    case ShortEdgeThresholdMode::HybridGlobalAndLocal:
        return setLengthThresholdHybrid();
    default:
        return false;
    }
}

double MergeShortEdge::computeStatisticBasedLengthThreshold()
{
    // 先在对数域做低分位统计，降低极长边对自动阈值的拉偏影响。
    //const std::vector<TopoDS_Edge> edges = collectAllEdges();
    if (allEdges_.empty())
        collectAllEdges();
    std::vector<double> logLengths;
    logLengths.reserve(allEdges_.size());

    for (const auto& edge : allEdges_) {
        const double length = GeometryUtils::ComputeEdgeLength(edge);
        if (length > Precision::Confusion() * 10.0) {
            logLengths.push_back(std::log(length));
        }
    }

    if (logLengths.size() < 5) {
        return 0.0;
    }

    std::sort(logLengths.begin(), logLengths.end());
    //取排序后的前10%
    const std::size_t quantileIndex = static_cast<std::size_t>(0.10 * static_cast<double>(logLengths.size() - 1));
    return std::exp(logLengths[quantileIndex]);
}

bool MergeShortEdge::satisfiesLocalLengthCriterion(const TopoDS_Edge& edge) const
{
    // 纯手工或纯 bbox 模式下保持原行为，只在混合模式引入局部相对长度约束。
    if (thresholdMode_ != ShortEdgeThresholdMode::HybridGlobalAndLocal) {
        return true;
    }

    if (vertexEdgeMap_.Extent() == 0) {
        return true;
    }

    const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
    //长度太短作为退化边处理，直接返回
    if (edgeLength <= Precision::Confusion()) {
        return true;
    }

    double minAdjacentLength = std::numeric_limits<double>::max();
    const TopoDS_Vertex vertices[] = {
        TopExp::FirstVertex(edge),
        TopExp::LastVertex(edge)
    };

    for (const auto& vertex : vertices) {
        if (vertex.IsNull() || !vertexEdgeMap_.Contains(vertex)) {
            continue;
        }

        const TopTools_ListOfShape& edgeList = vertexEdgeMap_.FindFromKey(vertex);
        for (TopTools_ListIteratorOfListOfShape it(edgeList); it.More(); it.Next()) {
            const TopoDS_Edge adjacent = TopoDS::Edge(it.Value());
            if (adjacent.IsNull() || adjacent.IsSame(edge)) {
                continue;
            }

            const double adjacentLength = GeometryUtils::ComputeEdgeLength(adjacent);
            if (adjacentLength > Precision::Confusion()) {
                minAdjacentLength = std::min(minAdjacentLength, adjacentLength);
            }
        }
    }

    if (minAdjacentLength == std::numeric_limits<double>::max()) {
        return true;
    }

    return edgeLength <= minAdjacentLength * localLengthRatio_;
}

bool MergeShortEdge::validateMergedEdgeGeometry(const TopoDS_Edge& edge) const
{
    // 这里只做轻量几何验收：至少保证端点仍能稳定投影回新边曲线。
    if (edge.IsNull()) {
        return false;
    }

    TopLoc_Location loc;
    Standard_Real first = 0.0;
    Standard_Real last = 0.0;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, loc, first, last);
    if (curve.IsNull()) {
        return false;
    }

    if (!loc.IsIdentity()) {
        curve = Handle(Geom_Curve)::DownCast(curve->Transformed(loc.Transformation()));
    }

    const TopoDS_Vertex vertices[] = {
        TopExp::FirstVertex(edge),
        TopExp::LastVertex(edge)
    };

    for (const auto& vertex : vertices) {
        if (vertex.IsNull()) {
            return false;
        }

        const gp_Pnt point = BRep_Tool::Pnt(vertex);
        GeomAPI_ProjectPointOnCurve projector(point, curve, first, last);
        if (projector.NbPoints() <= 0) {
            return false;
        }

        const double allowed = std::max({
            BRep_Tool::Tolerance(vertex),
            BRep_Tool::Tolerance(edge),
            Precision::Confusion()
            });
        if (projector.LowerDistance() > allowed) {
            return false;
        }
    }

    return true;
}

bool MergeShortEdge::hasValidPCurveOnFace(const TopoDS_Edge& edge, const TopoDS_Face& face) const
{
    return GeometryUtils::HasValidPCurveOnFace(edge, face);
}

bool MergeShortEdge::repairPCurvesForEdge(TopoDS_Edge& edge, const TopTools_ListOfShape& faces) const
{
    return GeometryUtils::RepairPCurvesForEdge(
        edge,
        faces,
        std::max(BRep_Tool::Tolerance(edge), Precision::Confusion()));
}

bool MergeShortEdge::repairAffectedPCurves(TopoDS_Shape& shape, std::vector<TopoDS_Edge>& newEdges) const
{
    if (shape.IsNull() || newEdges.empty()) {
        return true;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    bool allEdgesValid = true;
    for (auto& newEdge : newEdges) {
        if (newEdge.IsNull()) {
            allEdgesValid = false;
            continue;
        }

        if (!edgeFaceMap.Contains(newEdge)) {
            continue;
        }

        if (!repairPCurvesForEdge(newEdge, edgeFaceMap.FindFromKey(newEdge))) {
            allEdgesValid = false;
        }
    }

    return allEdgesValid;
}

bool MergeShortEdge::validateShape(const TopoDS_Shape& shape) const
{
    // merge 成功不等于模型有效，最终仍以 BRepCheck 的整体校验为准。
    if (shape.IsNull()) {
        return false;
    }

    BRepCheck_Analyzer analyzer(shape, Standard_True);
    return analyzer.IsValid();
}

Handle(Geom_Curve) MergeShortEdge::getUnderlyingCurve(const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return Handle(Geom_Curve)();
    }

    TopLoc_Location loc;
    Standard_Real first = 0.0;
    Standard_Real last = 0.0;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, loc, first, last);
    if (curve.IsNull()) {
        return curve;
    }

    if (!loc.IsIdentity()) {
        curve = Handle(Geom_Curve)::DownCast(curve->Transformed(loc.Transformation()));
    }

    while (curve->IsKind(STANDARD_TYPE(Geom_TrimmedCurve))) {
        curve = Handle(Geom_TrimmedCurve)::DownCast(curve)->BasisCurve();
    }

    return curve;
}

TopoDS_Edge MergeShortEdge::createMergedEdgeByCompositeBSpline(const MergeCandidate& candidate) const
{
    const TopoDS_Vertex v1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    const TopoDS_Vertex v2 = getOtherVertex(candidate.edge2, candidate.commonVertex);
    if (v1.IsNull() || v2.IsNull()) {
        return TopoDS_Edge();
    }

    try {
        auto makeTrimmedSegment = [&](const TopoDS_Edge& edge) -> Handle(Geom_TrimmedCurve) {
            TopLoc_Location loc;
            Standard_Real first = 0.0;
            Standard_Real last = 0.0;
            Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, loc, first, last);
            if (curve.IsNull()) {
                return Handle(Geom_TrimmedCurve)();
            }

            if (!loc.IsIdentity()) {
                curve = Handle(Geom_Curve)::DownCast(curve->Transformed(loc.Transformation()));
            }
            if (curve.IsNull()) {
                return Handle(Geom_TrimmedCurve)();
            }

            return new Geom_TrimmedCurve(curve, first, last);
        };

        Handle(Geom_TrimmedCurve) edge1Seg = makeTrimmedSegment(candidate.edge1);
        Handle(Geom_TrimmedCurve) edge2Seg = makeTrimmedSegment(candidate.edge2);
        if (edge1Seg.IsNull() || edge2Seg.IsNull()) {
            return TopoDS_Edge();
        }

        auto orientSegment = [&](const Handle(Geom_TrimmedCurve)& segment,
                                 const TopoDS_Edge& edge,
                                 const TopoDS_Vertex& startVertex,
                                 const TopoDS_Vertex& endVertex) -> bool {
            if (segment.IsNull() || startVertex.IsNull() || endVertex.IsNull()) {
                return false;
            }

            const TopoDS_Vertex edgeFirst = TopExp::FirstVertex(edge);
            const TopoDS_Vertex edgeLast = TopExp::LastVertex(edge);
            if (edgeFirst.IsSame(startVertex) && edgeLast.IsSame(endVertex)) {
                return true;
            }
            if (edgeFirst.IsSame(endVertex) && edgeLast.IsSame(startVertex)) {
                segment->Reverse();
                return true;
            }

            return false;
        };

        Handle(Geom_TrimmedCurve) firstSeg;
        Handle(Geom_TrimmedCurve) secondSeg;
        TopoDS_Vertex startVertex;
        TopoDS_Vertex endVertex;

        // 复合顺序以 edge2 的方向为主导：
        // 若 edge2 从公共点出发，则整体方向为 other(edge1) -> common -> other(edge2)；
        // 若 edge2 终止于公共点，则整体方向为 other(edge2) -> common -> other(edge1)。
        if (candidate.commonVertex.IsSame(TopExp::FirstVertex(candidate.edge2))) {
            firstSeg = edge1Seg;
            secondSeg = edge2Seg;
            startVertex = v1;
            endVertex = v2;

            if (!orientSegment(firstSeg, candidate.edge1, v1, candidate.commonVertex) ||
                !orientSegment(secondSeg, candidate.edge2, candidate.commonVertex, v2)) {
                return TopoDS_Edge();
            }
        }
        else if (candidate.commonVertex.IsSame(TopExp::LastVertex(candidate.edge2))) {
            firstSeg = edge2Seg;
            secondSeg = edge1Seg;
            startVertex = v2;
            endVertex = v1;

            if (!orientSegment(firstSeg, candidate.edge2, v2, candidate.commonVertex) ||
                !orientSegment(secondSeg, candidate.edge1, candidate.commonVertex, v1)) {
                return TopoDS_Edge();
            }
        }
        else {
            return TopoDS_Edge();
        }

        GeomConvert_CompCurveToBSplineCurve comp(firstSeg, Convert_QuasiAngular);
        const Standard_Real tol = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(candidate.commonVertex),
            BRep_Tool::Tolerance(candidate.edge1),
            BRep_Tool::Tolerance(candidate.edge2)
            });
        if (!comp.Add(secondSeg, tol, Standard_True, Standard_True, 2)) {
            return TopoDS_Edge();
        }

        Handle(Geom_BSplineCurve) mergedCurve = comp.BSplineCurve();
        if (mergedCurve.IsNull()) {
            return TopoDS_Edge();
        }

        BRepBuilderAPI_MakeEdge edgeMaker(mergedCurve, startVertex, endVertex);
        if (edgeMaker.IsDone()) {
            TopoDS_Edge newEdge = edgeMaker.Edge();
            newEdge.Orientation(candidate.edge2.Orientation());
            return newEdge;
        }
    }
    catch (Standard_Failure&) {
        std::cout << "MergedEdge Failed!" << std::endl;
    }

    return TopoDS_Edge();
}

std::vector<TopoDS_Edge> MergeShortEdge::getShortEdges()
{
    if (shortEdges_.empty()) {
        detectShortEdges();
    }
    return shortEdges_;
}

void MergeShortEdge::collectAllEdges()
{
    allEdges_.clear();
    if (inputShape_.IsNull()) {
        return;
    }

    TopTools_MapOfShape addedEdges;

    for (TopExp_Explorer exp(inputShape_, TopAbs_EDGE); exp.More(); exp.Next()) {
        const TopoDS_Edge currentEdge = TopoDS::Edge(exp.Current());
        if (addedEdges.Add(currentEdge)) {
            allEdges_.push_back(currentEdge);
        }
    }
}

bool MergeShortEdge::isShortEdge(const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return false;
    }

    const double length = GeometryUtils::ComputeEdgeLength(edge);
    return (length <= lengthThreshold_ + Precision::Confusion()) &&
        satisfiesLocalLengthCriterion(edge);
}

bool MergeShortEdge::isHighRiskEdge(const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return true;
    }

    if (BRep_Tool::Degenerated(edge)) {
        return true;
    }

    if (BRep_Tool::IsClosed(edge)) {
        return true;
    }

    TopLoc_Location loc;
    Standard_Real first = 0.0;
    Standard_Real last = 0.0;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, loc, first, last);
    if (curve.IsNull()) {
        return true;
    }

    const auto hasSeamOnFaces = [&](const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap) {
        if (!edgeFaceMap.Contains(edge)) {
            return false;
        }

        const TopTools_ListOfShape& faces = edgeFaceMap.FindFromKey(edge);
        for (TopTools_ListIteratorOfListOfShape it(faces); it.More(); it.Next()) {
            const TopoDS_Face face = TopoDS::Face(it.Value());
            if (!face.IsNull() && BRep_Tool::IsClosed(edge, face)) {
                return true;
            }
        }

        return false;
    };

    if (hasSeamOnFaces(edgeFaceMap_)) {
        return true;
    }

    return false;
}

bool MergeShortEdge::findBestAdjacentEdge(const TopoDS_Edge& shortEdge, TopoDS_Vertex& bestCommonVertex, TopoDS_Edge& bestAdjacentEdge) const
{
    if (isHighRiskEdge(shortEdge)) {
        return false;
    }

    const TopoDS_Vertex v1 = TopExp::FirstVertex(shortEdge);
    const TopoDS_Vertex v2 = TopExp::LastVertex(shortEdge);

    if (v1.IsNull() || v2.IsNull()) {
        return false;
    }

    double maxAngle = -1.0;
    bool found = false;

    auto evaluateVertex = [&](const TopoDS_Vertex& v) {
        if (!vertexEdgeMap_.Contains(v)) {
            return;
        }

        const TopTools_ListOfShape& edgeList = vertexEdgeMap_.FindFromKey(v);
        if (edgeList.Extent() != 2) {
            return;
        }

        for (TopTools_ListIteratorOfListOfShape it(edgeList); it.More(); it.Next()) {
            const TopoDS_Edge adjEdge = TopoDS::Edge(it.Value());
            if (adjEdge.IsSame(shortEdge)) {
                continue;
            }

            if (isHighRiskEdge(adjEdge)) {
                continue;
            }

            if (isShortEdge(adjEdge)) {
                continue;
            }

            const double angle = GeometryUtils::ComputeIntersectionAngle(shortEdge, adjEdge, v);
            const double angleDeg = angle * 180.0 / M_PI;
            if (angleDeg >= angleThresholdDeg_ && angleDeg > maxAngle) {
                maxAngle = angleDeg;
                bestCommonVertex = v;
                bestAdjacentEdge = adjEdge;
                found = true;
            }
        }
    };

    evaluateVertex(v1);
    evaluateVertex(v2);
    return found;
}

void MergeShortEdge::detectShortEdges()
{
    shortEdges_.clear();

    if (inputShape_.IsNull()) {
        return;
    }

    if (vertexEdgeMap_.Extent() == 0 || edgeFaceMap_.Extent() == 0) {
        buildTopologyMap();
    }

    //const std::vector<TopoDS_Edge> allEdges = collectAllEdges();
    if (allEdges_.empty())
        collectAllEdges();

    for (const auto& edge : allEdges_) {
        if (isHighRiskEdge(edge)) {
            continue;
        }

        if (isShortEdge(edge)) {
            shortEdges_.push_back(edge);
        }
    }
}

TopoDS_Shape MergeShortEdge::preprocessWithShapeFixWireframe(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return shape;
    }

    ShapeFix_Wireframe wireframeFix(shape);
    const double preFixPrecision =
        (lengthThreshold_ > 0.0)
        ? std::max(Precision::Confusion(), lengthThreshold_ * 0.5)
        : Precision::Confusion();

    wireframeFix.SetPrecision(preFixPrecision);
    wireframeFix.ModeDropSmallEdges() = Standard_False;
    
    wireframeFix.SetLimitAngle(angleThresholdDeg_ * M_PI / 180.0);
    wireframeFix.FixSmallEdges();

    const TopoDS_Shape fixedShape = wireframeFix.Shape();
    return fixedShape.IsNull() ? shape : fixedShape;
}

std::vector<MergeCandidate> MergeShortEdge::findMergeCandidates()
{
    mergeCandidates_.clear();
    if (inputShape_.IsNull()) {
        return mergeCandidates_;
    }

    buildTopologyMap();
    if (shortEdges_.empty()) {
        detectShortEdges();
    }

    TopTools_MapOfShape processedEdges;

    for (const auto& shortEdge : shortEdges_) {
        if (processedEdges.Contains(shortEdge)) {
            continue;
        }

        TopoDS_Vertex commonVertex;
        TopoDS_Edge adjacentEdge;

        if (!findBestAdjacentEdge(shortEdge, commonVertex, adjacentEdge)) {
            continue;
        }

        if (processedEdges.Contains(adjacentEdge)) {
            continue;
        }

        MergeCandidate candidate;
        candidate.edge1 = shortEdge;
        candidate.edge2 = adjacentEdge;
        candidate.commonVertex = commonVertex;
        candidate.length1 = GeometryUtils::ComputeEdgeLength(shortEdge);
        candidate.length2 = GeometryUtils::ComputeEdgeLength(adjacentEdge);
        candidate.angle = GeometryUtils::ComputeIntersectionAngle(shortEdge, adjacentEdge, commonVertex);

        mergeCandidates_.push_back(candidate);

        processedEdges.Add(shortEdge);
        processedEdges.Add(adjacentEdge);
    }

    return mergeCandidates_;
}

int MergeShortEdge::performMerge()
{
    if (inputShape_.IsNull()) {
        return 0;
    }

    //resultShape_ = inputShape_;
    //resultShape_ = preprocessWithShapeFixWireframe(inputShape_);
    //inputShape_ = resultShape_;
    int mergeCount = 0;

    refreshLengthThreshold();
    allEdges_.clear();
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    edgeFaceMap_.Clear();

    findMergeCandidates();
    if (mergeCandidates_.empty()) {
        return 0;
    }

    const std::vector<std::vector<MergeCandidate>> candidateBatches = buildCandidateBatches(mergeCandidates_);
    TopoDS_Shape currentShape = inputShape_;

    for (const auto& batch : candidateBatches) {
        TopoDS_Shape committedShape;
        int committedCount = 0;
        if (!tryApplyBatch(currentShape, batch, committedShape, committedCount)) {
            continue;
        }

        currentShape = committedShape;
        mergeCount += committedCount;
    }

    resultShape_ = currentShape;
    inputShape_ = resultShape_;
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    edgeFaceMap_.Clear();
    refreshLengthThreshold();

    return mergeCount;
}

std::vector<std::vector<MergeCandidate>> MergeShortEdge::buildCandidateBatches(const std::vector<MergeCandidate>& candidates) const
{
    std::vector<std::vector<MergeCandidate>> batches;
    // 第 i 个 batch 里已经占用过的 edge 集合
    std::vector<TopTools_MapOfShape> batchEdges;
    // 第 i 个 batch 里已经占用过的 vertex 集合
    std::vector<TopTools_MapOfShape> batchVertices;

    for (const auto& candidate : candidates) {
        // 取出每个candidate会影响到的边和顶点
        const TopoDS_Vertex otherVertex1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
        const TopoDS_Vertex otherVertex2 = getOtherVertex(candidate.edge2, candidate.commonVertex);

        const TopoDS_Shape candidateEdges[] = { candidate.edge1, candidate.edge2 };
        const TopoDS_Shape candidateVertices[] = { candidate.commonVertex, otherVertex1, otherVertex2 };

        bool placed = false;
        for (std::size_t i = 0; i < batches.size(); ++i) {
            bool hasConflict = false;

            for (const auto& edge : candidateEdges) {
                if (!edge.IsNull() && batchEdges[i].Contains(edge)) {
                    hasConflict = true;
                    break;
                }
            }

            if (!hasConflict) {
                for (const auto& vertex : candidateVertices) {
                    if (!vertex.IsNull() && batchVertices[i].Contains(vertex)) {
                        hasConflict = true;
                        break;
                    }
                }
            }

            if (hasConflict) {
                continue;
            }

            batches[i].push_back(candidate);
            for (const auto& edge : candidateEdges) {
                if (!edge.IsNull()) {
                    batchEdges[i].Add(edge);
                }
            }
            for (const auto& vertex : candidateVertices) {
                if (!vertex.IsNull()) {
                    batchVertices[i].Add(vertex);
                }
            }
            placed = true;
            break;
        }

        if (placed) {
            continue;
        }

        batches.emplace_back(1, candidate);
        batchEdges.emplace_back();
        batchVertices.emplace_back();
        for (const auto& edge : candidateEdges) {
            if (!edge.IsNull()) {
                batchEdges.back().Add(edge);
            }
        }
        for (const auto& vertex : candidateVertices) {
            if (!vertex.IsNull()) {
                batchVertices.back().Add(vertex);
            }
        }
    }

    return batches;
}

bool MergeShortEdge::tryApplyBatch(const TopoDS_Shape& baseShape,
                                   const std::vector<MergeCandidate>& batch,
                                   TopoDS_Shape& committedShape,
                                   int& committedCount) const
{
    committedShape = TopoDS_Shape();
    committedCount = 0;

    if (baseShape.IsNull() || batch.empty()) {
        return false;
    }

    BRep_Builder builder;
    BRepTools_ReShape batchReshaper;
    std::vector<TopoDS_Edge> batchNewEdges;
    TopTools_MapOfShape batchRemovedEdges;

    for (const auto& candidate : batch) {
        if (batchRemovedEdges.Contains(candidate.edge1) || batchRemovedEdges.Contains(candidate.edge2)) {
            continue;
        }

        TopoDS_Edge newEdge = createMergedEdge(candidate);
        if (newEdge.IsNull() || !validateMergedEdgeGeometry(newEdge)) {
            continue;
        }

        const double safeTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(candidate.commonVertex),
            candidate.length1 * 0.5
            });
        builder.UpdateEdge(newEdge, safeTolerance);

        batchReshaper.Replace(candidate.edge2, newEdge);
        batchReshaper.Remove(candidate.edge1);

        batchNewEdges.push_back(newEdge);
        batchRemovedEdges.Add(candidate.edge1);
        batchRemovedEdges.Add(candidate.edge2);
        committedCount++;
    }

    if (committedCount == 0) {
        return false;
    }

    TopoDS_Shape trialShape = batchReshaper.Apply(baseShape);
    if (trialShape.IsNull()) {
        committedCount = 0;
        return false;
    }

    if (!repairAffectedPCurves(trialShape, batchNewEdges)) {
        committedCount = 0;
        return false;
    }

    ShapeFix_Shape fixer(trialShape);
    const double fixPrecision = std::max(lengthThreshold_ * 2.0, Precision::Confusion());
    fixer.SetPrecision(fixPrecision);
    fixer.SetMaxTolerance(std::max(lengthThreshold_ * 5.0, fixPrecision));
    fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
    fixer.FixFaceTool()->FixMissingSeamMode() = 1;
    fixer.Perform();
    trialShape = fixer.Shape();

    if (!validateShape(trialShape)) {
        committedCount = 0;
        return false;
    }

    committedShape = trialShape;
    return true;
}

TopoDS_Edge MergeShortEdge::createMergedEdge(const MergeCandidate& candidate) const
{
    const TopoDS_Vertex v1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    const TopoDS_Vertex v2 = getOtherVertex(candidate.edge2, candidate.commonVertex);

    if (v1.IsNull() || v2.IsNull()) {
        return TopoDS_Edge();
    }

    const gp_Pnt p1 = BRep_Tool::Pnt(v1);
    const gp_Pnt p2 = BRep_Tool::Pnt(v2);
    const gp_Pnt pCommon = BRep_Tool::Pnt(candidate.commonVertex);

    Handle(Geom_Curve) baseCurve = getUnderlyingCurve(candidate.edge2);
    if (!baseCurve.IsNull()) {
        TopoDS_Edge compositeSplineEdge = createMergedEdgeByCompositeBSpline(candidate);
        if (!compositeSplineEdge.IsNull()) {
            return compositeSplineEdge;
        }

        try {
            if (baseCurve->IsKind(STANDARD_TYPE(Geom_Line))) {
                BRepBuilderAPI_MakeEdge edgeMaker(v1, v2);
                if (edgeMaker.IsDone()) {
                    TopoDS_Edge newEdge = edgeMaker.Edge();
                    newEdge.Orientation(candidate.edge2.Orientation());
                    return newEdge;
                }
            }
            else if (baseCurve->IsKind(STANDARD_TYPE(Geom_Circle))) {
                GeomAPI_ProjectPointOnCurve projector(p1, baseCurve);
                gp_Pnt projectedP1 = p1;
                if (projector.NbPoints() > 0) {
                    projectedP1 = projector.NearestPoint();
                }

                GC_MakeArcOfCircle arcMaker(projectedP1, pCommon, p2);
                if (arcMaker.IsDone()) {
                    BRepBuilderAPI_MakeEdge edgeMaker(arcMaker.Value(), v1, v2);
                    if (edgeMaker.IsDone()) {
                        TopoDS_Edge newEdge = edgeMaker.Edge();
                        newEdge.Orientation(candidate.edge2.Orientation());
                        return newEdge;
                    }
                }
            }
            else if (baseCurve->IsKind(STANDARD_TYPE(Geom_Ellipse)) ||
                baseCurve->IsKind(STANDARD_TYPE(Geom_Parabola)) ||
                baseCurve->IsKind(STANDARD_TYPE(Geom_Hyperbola))) {
                BRepBuilderAPI_MakeEdge edgeMaker(baseCurve, v1, v2);
                if (edgeMaker.IsDone()) {
                    TopoDS_Edge newEdge = edgeMaker.Edge();
                    newEdge.Orientation(candidate.edge2.Orientation());
                    return newEdge;
                }
            }
            else if (baseCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)) ||
                baseCurve->IsKind(STANDARD_TYPE(Geom_BezierCurve))) {
                Handle(Geom_BoundedCurve) boundedCurve = Handle(Geom_BoundedCurve)::DownCast(baseCurve->Copy());

                if (!boundedCurve.IsNull()) {
                    const Standard_Real uCommon = BRep_Tool::Parameter(candidate.commonVertex, candidate.edge2);
                    const Standard_Real bcFirst = boundedCurve->FirstParameter();
                    const Standard_Real bcLast = boundedCurve->LastParameter();
                    const Standard_Boolean extendAfter = (Abs(uCommon - bcLast) < Abs(uCommon - bcFirst));

                    GeomAPI_ProjectPointOnCurve projP1(p1, boundedCurve);
                    const Standard_Real tolV1 = BRep_Tool::Tolerance(v1);
                    bool needsExtension = true;

                    if (projP1.NbPoints() > 0) {
                        const Standard_Real uProj = projP1.LowerDistanceParameter();
                        if (uProj >= bcFirst - Precision::Confusion() &&
                            uProj <= bcLast + Precision::Confusion() &&
                            projP1.LowerDistance() <= tolV1) {
                            needsExtension = false;
                        }
                    }

                    if (needsExtension) {
                        bool extended = false;
                        try {
                            GeomLib::ExtendCurveToPoint(boundedCurve, p1, 1, extendAfter);
                            extended = true;
                        }
                        catch (Standard_Failure&) {
                            try {
                                GeomLib::ExtendCurveToPoint(boundedCurve, p1, 0, extendAfter);
                                extended = true;
                            }
                            catch (Standard_Failure&) {
                            }
                        }

                        if (!extended) {
                            Standard_Failure::Raise("B-Spline extension failed");
                        }
                    }

                    BRepBuilderAPI_MakeEdge edgeMaker(boundedCurve, v1, v2);
                    if (edgeMaker.IsDone()) {
                        TopoDS_Edge newEdge = edgeMaker.Edge();
                        newEdge.Orientation(candidate.edge2.Orientation());
                        return newEdge;
                    }
                }
            }
        }
        catch (Standard_Failure&) {
        }
    }

    try {
        BRepBuilderAPI_MakeEdge fallbackMaker(v1, v2);
        if (fallbackMaker.IsDone()) {
            TopoDS_Edge newEdge = fallbackMaker.Edge();
            newEdge.Orientation(candidate.edge2.Orientation());
            return newEdge;
        }
    }
    catch (Standard_Failure&) {
        return TopoDS_Edge();
    }

    return TopoDS_Edge();
}

TopoDS_Vertex MergeShortEdge::getOtherVertex(const TopoDS_Edge& edge, const TopoDS_Vertex& commonVertex) const
{
    if (edge.IsNull() || commonVertex.IsNull()) {
        return TopoDS_Vertex();
    }

    const TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
    const TopoDS_Vertex v2 = TopExp::LastVertex(edge);
    if (v1.IsSame(commonVertex)) {
        return v2;
    }
    if (v2.IsSame(commonVertex)) {
        return v1;
    }
    return TopoDS_Vertex();
}

bool MergeShortEdge::setLengthThresholdByBoundingBox(double ratio)
{
    if (inputShape_.IsNull() || ratio <= 0.0) {
        return false;
    }

    if (ratio > 0.01) {
        ratio = 0.01;
    }

    const double diagonal = GeometryUtils::ComputeBoundingBoxDiagonal(inputShape_);
    if (diagonal <= 0.0) {
        return false;
    }

    bboxRatio_ = ratio;
    lengthThreshold_ = diagonal * ratio;
    thresholdMode_ = ShortEdgeThresholdMode::BoundingBoxRatio;
    return true;
}

bool MergeShortEdge::setLengthThresholdHybrid()
{
    if (inputShape_.IsNull()) {
        return false;
    }

    const double statThreshold = computeStatisticBasedLengthThreshold();
    double bboxThreshold = 0.0;
    if (bboxRatio_ > 0.0) {
        const double diagonal = GeometryUtils::ComputeBoundingBoxDiagonal(inputShape_);
        if (diagonal > 0.0) {
            bboxThreshold = diagonal * std::min(bboxRatio_, 0.01);
        }
    }

    // 统计阈值负责反映真实边长分布，bbox 阈值只作为保守上限和兜底。
    double resolvedThreshold = statThreshold;
    if (resolvedThreshold <= 0.0) {
        resolvedThreshold = bboxThreshold;
    }
    else if (bboxThreshold > 0.0) {
        resolvedThreshold = std::min(resolvedThreshold, bboxThreshold);
    }

    if (resolvedThreshold <= 0.0) {
        return false;
    }

    lengthThreshold_ = std::max(resolvedThreshold, Precision::Confusion() * 10.0);
    thresholdMode_ = ShortEdgeThresholdMode::HybridGlobalAndLocal;
    return true;
}
