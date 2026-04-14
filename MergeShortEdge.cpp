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
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepTools_ReShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <gp_Pnt.hxx>
#include <TopTools_MapOfShape.hxx>
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
    inputShapes_.clear();
    if (!shape.IsNull()) {
        inputShapes_.push_back(shape);
    }

    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    resultShape_.Nullify();

    refreshLengthThreshold();
}

void MergeShortEdge::setInputShapes(const std::vector<TopoDS_Shape>& shapes)
{
    inputShapes_ = shapes;
    std::erase_if(inputShapes_, [](const TopoDS_Shape& s) { return s.IsNull(); });

    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    resultShape_.Nullify();

    refreshLengthThreshold();
}

void MergeShortEdge::buildTopologyMap()
{
    vertexEdgeMap_.Clear();

    TopoDS_Compound compound;
    BRep_Builder builder;
    builder.MakeCompound(compound);
    for (const auto& shape : inputShapes_) {
        builder.Add(compound, shape);
    }

    TopologyUtils::BuildVertexEdgeMap(compound, vertexEdgeMap_);
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

double MergeShortEdge::computeStatisticBasedLengthThreshold() const
{
    // 先在对数域做低分位统计，降低极长边对自动阈值的拉偏影响。
    const std::vector<TopoDS_Edge> edges = collectAllEdges();
    std::vector<double> logLengths;
    logLengths.reserve(edges.size());

    for (const auto& edge : edges) {
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
    const std::size_t quantileIndex = static_cast<std::size_t>(0.1 * static_cast<double>(logLengths.size() - 1));
    return std::exp(logLengths[quantileIndex]) * 0.5;
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
    if (edge.IsNull() || face.IsNull()) {
        return false;
    }

    Standard_Real first = 0.0;
    Standard_Real last = 0.0;
    Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
    return !pcurve.IsNull();
}

bool MergeShortEdge::repairPCurvesForEdge(TopoDS_Edge& edge, const TopTools_ListOfShape& faces) const
{
    if (edge.IsNull()) {
        return false;
    }

    ShapeFix_Edge edgeFix;
    bool allFacesValid = true;

    for (TopTools_ListIteratorOfListOfShape it(faces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        if (!hasValidPCurveOnFace(edge, face)) {
            const Standard_Boolean isSeam = BRep_Tool::IsClosed(edge, face);
            edgeFix.FixAddPCurve(edge, face, isSeam, std::max(BRep_Tool::Tolerance(edge), Precision::Confusion()));
        }

        edgeFix.FixReversed2d(edge, face);
        edgeFix.FixSameParameter(edge, face, std::max(BRep_Tool::Tolerance(edge), Precision::Confusion()));
        edgeFix.FixVertexTolerance(edge, face);

        if (!hasValidPCurveOnFace(edge, face)) {
            allFacesValid = false;
        }
    }

    return allFacesValid;
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

TopoDS_Edge MergeShortEdge::createMergedEdgeFromSharedBasisCurve(
    const MergeCandidate& candidate,
    const Handle(Geom_Curve)& sharedCurve) const
{
    if (sharedCurve.IsNull()) {
        return TopoDS_Edge();
    }

    const TopoDS_Vertex v1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    const TopoDS_Vertex v2 = getOtherVertex(candidate.edge2, candidate.commonVertex);
    if (v1.IsNull() || v2.IsNull()) {
        return TopoDS_Edge();
    }

    try {
        BRepBuilderAPI_MakeEdge edgeMaker(sharedCurve, v1, v2);
        if (edgeMaker.IsDone()) {
            TopoDS_Edge newEdge = edgeMaker.Edge();
            newEdge.Orientation(candidate.edge2.Orientation());
            return newEdge;
        }
    }
    catch (Standard_Failure&) {
    }

    return TopoDS_Edge();
}

TopoDS_Edge MergeShortEdge::createMergedEdgeByCompositeBSpline(const MergeCandidate& candidate) const
{
    const TopoDS_Vertex v1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    const TopoDS_Vertex v2 = getOtherVertex(candidate.edge2, candidate.commonVertex);
    if (v1.IsNull() || v2.IsNull()) {
        return TopoDS_Edge();
    }

    try {
        BRepAdaptor_Curve adaptor1(candidate.edge1);
        BRepAdaptor_Curve adaptor2(candidate.edge2);
        Handle(Geom_BSplineCurve) bs1 = adaptor1.BSpline();
        Handle(Geom_BSplineCurve) bs2 = adaptor2.BSpline();
        if (bs1.IsNull() || bs2.IsNull()) {
            return TopoDS_Edge();
        }

        Handle(Geom_BSplineCurve) seg1 = Handle(Geom_BSplineCurve)::DownCast(bs1->Copy());
        Handle(Geom_BSplineCurve) seg2 = Handle(Geom_BSplineCurve)::DownCast(bs2->Copy());
        if (seg1.IsNull() || seg2.IsNull()) {
            return TopoDS_Edge();
        }
		//确保两段底层曲线朝向一致
        if (candidate.commonVertex.IsSame(TopExp::FirstVertex(candidate.edge1))) {
            seg1->Reverse();
        }
        if (candidate.commonVertex.IsSame(TopExp::LastVertex(candidate.edge2))) {
            seg2->Reverse();
        }
        
        GeomConvert_CompCurveToBSplineCurve comp(seg1);
        const Standard_Real tol = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(candidate.commonVertex),
            BRep_Tool::Tolerance(candidate.edge1),
            BRep_Tool::Tolerance(candidate.edge2)
            });
        if (!comp.Add(seg2, tol, Standard_True, Standard_True, 0)) {
            return TopoDS_Edge();
        }

        Handle(Geom_BSplineCurve) mergedCurve = comp.BSplineCurve();
        if (mergedCurve.IsNull()) {
            return TopoDS_Edge();
        }

        BRepBuilderAPI_MakeEdge edgeMaker(mergedCurve, v1, v2);
        if (edgeMaker.IsDone()) {
            TopoDS_Edge newEdge = edgeMaker.Edge();
            newEdge.Orientation(candidate.edge2.Orientation());
            return newEdge;
        }
    }
    catch (Standard_Failure&) {
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

std::vector<TopoDS_Edge> MergeShortEdge::collectAllEdges() const
{
    std::vector<TopoDS_Edge> edges;
    TopTools_MapOfShape addedEdges;

    for (const auto& shape : inputShapes_) {
        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            const TopoDS_Edge currentEdge = TopoDS::Edge(exp.Current());
            if (addedEdges.Add(currentEdge)) {
                edges.push_back(currentEdge);
            }
        }
    }
    return edges;
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

bool MergeShortEdge::findBestAdjacentEdge(const TopoDS_Edge& shortEdge, TopoDS_Vertex& bestCommonVertex, TopoDS_Edge& bestAdjacentEdge) const
{
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

    if (inputShapes_.empty()) {
        return;
    }

    if (thresholdMode_ == ShortEdgeThresholdMode::HybridGlobalAndLocal &&
        vertexEdgeMap_.Extent() == 0) {
        // 混合模式依赖局部邻边长度，因此短边检测前要先具备顶点-边映射。
        buildTopologyMap();
    }

    const std::vector<TopoDS_Edge> allEdges = collectAllEdges();
    for (const auto& edge : allEdges) {
        if (isShortEdge(edge)) {
            shortEdges_.push_back(edge);
        }
    }
}

TopoDS_Shape MergeShortEdge::executeUnifySameDomain(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) {
        return shape;
    }

    ShapeUpgrade_UnifySameDomain unifier(shape, Standard_True, Standard_True, Standard_False);
    unifier.Build();
    return unifier.Shape();
}

std::vector<MergeCandidate> MergeShortEdge::findMergeCandidates()
{
    mergeCandidates_.clear();
    if (inputShapes_.empty()) {
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
    if (inputShapes_.empty()) {
        return 0;
    }

    BRep_Builder builder;
    TopoDS_Compound initialCompound;
    builder.MakeCompound(initialCompound);
    for (const auto& shape : inputShapes_) {
        builder.Add(initialCompound, shape);
    }

    resultShape_ = executeUnifySameDomain(initialCompound);
    inputShapes_.clear();
    inputShapes_.push_back(resultShape_);
    refreshLengthThreshold();

    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();

    findMergeCandidates();
    if (mergeCandidates_.empty()) {
        return 0;
    }

    const TopoDS_Shape baseShape = resultShape_;
    BRepTools_ReShape reshaper;
    int mergeCount = 0;
    std::vector<TopoDS_Edge> newEdges;
    TopTools_MapOfShape removedEdges;

    for (const auto& candidate : mergeCandidates_) {
        if (removedEdges.Contains(candidate.edge1) || removedEdges.Contains(candidate.edge2)) {
            continue;
        }

        TopoDS_Edge newEdge = createMergedEdge(candidate);
        if (newEdge.IsNull() || !validateMergedEdgeGeometry(newEdge)) {
            // 局部重建不可靠时直接跳过，不把问题留给后面的 ShapeFix 去兜底。
            continue;
        }

        const double safeTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(candidate.commonVertex),
            candidate.length1 * 0.5
            });
        builder.UpdateEdge(newEdge, safeTolerance);

        reshaper.Replace(candidate.edge2, newEdge);
        reshaper.Remove(candidate.edge1);

        newEdges.push_back(newEdge);
        removedEdges.Add(candidate.edge1);
        removedEdges.Add(candidate.edge2);
        mergeCount++;
    }

    if (mergeCount == 0) {
        resultShape_ = baseShape;
        return 0;
    }

    resultShape_ = reshaper.Apply(resultShape_);

    if (!repairAffectedPCurves(resultShape_, newEdges)) {
        resultShape_ = baseShape;
        inputShapes_.clear();
        inputShapes_.push_back(resultShape_);
        shortEdges_.clear();
        mergeCandidates_.clear();
        vertexEdgeMap_.Clear();
        refreshLengthThreshold();
        return 0;
    }

    ShapeFix_Shape fixer(resultShape_);
    const double fixPrecision = std::max(lengthThreshold_ * 2.0, Precision::Confusion());
    fixer.SetPrecision(fixPrecision);
    fixer.SetMaxTolerance(std::max(lengthThreshold_ * 5.0, fixPrecision));
    fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
    fixer.FixFaceTool()->FixMissingSeamMode() = 1;
    fixer.Perform();
    resultShape_ = fixer.Shape();

    if (!validateShape(resultShape_)) {
        // 整体验证失败时回退到 merge 前结果，避免输出一个“修过但更坏”的模型。
        resultShape_ = baseShape;
        inputShapes_.clear();
        inputShapes_.push_back(resultShape_);
        shortEdges_.clear();
        mergeCandidates_.clear();
        vertexEdgeMap_.Clear();
        refreshLengthThreshold();
        return 0;
    }

    inputShapes_.clear();
    inputShapes_.push_back(resultShape_);
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    refreshLengthThreshold();

    return mergeCount;
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
        const Handle(Geom_Curve) shortBasisCurve = getUnderlyingCurve(candidate.edge1);
        if (!shortBasisCurve.IsNull() && shortBasisCurve == baseCurve) {
            TopoDS_Edge sharedBasisEdge = createMergedEdgeFromSharedBasisCurve(candidate, baseCurve);
            if (!sharedBasisEdge.IsNull()) {
                return sharedBasisEdge;
            }
        }

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
    if (inputShapes_.empty() || ratio <= 0.0) {
        return false;
    }

    if (ratio > 0.01) {
        ratio = 0.01;
    }

    const double diagonal = GeometryUtils::ComputeBoundingBoxDiagonal(inputShapes_);
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
    if (inputShapes_.empty()) {
        return false;
    }

    const double statThreshold = computeStatisticBasedLengthThreshold();
    double bboxThreshold = 0.0;
    if (bboxRatio_ > 0.0) {
        const double diagonal = GeometryUtils::ComputeBoundingBoxDiagonal(inputShapes_);
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
