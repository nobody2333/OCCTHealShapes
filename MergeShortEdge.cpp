#include <algorithm>
#include <cmath>

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
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <gp_Pnt.hxx>
#include <TopTools_MapOfShape.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Parabola.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GC_MakeArcOfCircle.hxx>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MergeShortEdge::MergeShortEdge(double bboxRatio, double angleThresholdDeg)
    : lengthThreshold_(0.0)
    , angleThresholdDeg_(angleThresholdDeg)
    , bboxRatio_(bboxRatio)
{
}

void MergeShortEdge::setInputShape(const TopoDS_Shape& shape)
{
    inputShapes_.clear();
    if (!shape.IsNull()) {
        inputShapes_.push_back(shape);
    }
    // 重置结果
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    
    // 根据包围盒比例自动计算长度阈值
    if (bboxRatio_ > 0.0) {
        setLengthThresholdByBoundingBox(bboxRatio_);
    }
}

void MergeShortEdge::setInputShapes(const std::vector<TopoDS_Shape>& shapes)
{
    inputShapes_ = shapes;
    // 过滤掉空形状
    std::erase_if(inputShapes_, [](const TopoDS_Shape& s) {return s.IsNull(); });
    
    // 重置结果
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
    
    // 根据包围盒比例自动计算长度阈值
    if (bboxRatio_ > 0.0) {
        setLengthThresholdByBoundingBox(bboxRatio_);
    }
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

std::vector<TopoDS_Edge> MergeShortEdge::getShortEdges() {
    if(shortEdges_.empty()) {
        detectShortEdges();
	}
	return shortEdges_;
}

std::vector<TopoDS_Edge> MergeShortEdge::collectAllEdges() const
{
    std::vector<TopoDS_Edge> edges;
    TopTools_MapOfShape addedEdges; // 防止多重面共用边导致重复

    for (const auto& shape : inputShapes_) {
        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            TopoDS_Edge currentEdge = TopoDS::Edge(exp.Current());
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

    double length = GeometryUtils::ComputeEdgeLength(edge);
    
    // 使用 Precision::Confusion() 作为数值容差
    return (length <= lengthThreshold_ + Precision::Confusion());
}

bool MergeShortEdge::findBestAdjacentEdge(const TopoDS_Edge& shortEdge, TopoDS_Vertex& bestCommonVertex, TopoDS_Edge& bestAdjacentEdge) const
{
    TopoDS_Vertex v1 = TopExp::FirstVertex(shortEdge);
    TopoDS_Vertex v2 = TopExp::LastVertex(shortEdge);

    if (v1.IsNull() || v2.IsNull()) 
        return false;

    double maxAngle = -1.0;
    bool found = false;

    auto evaluateVertex = [&](const TopoDS_Vertex& v) {
        if (!vertexEdgeMap_.Contains(v)) 
            return;

        const TopTools_ListOfShape& edgeList = vertexEdgeMap_.FindFromKey(v);
        if (edgeList.Extent() != 2) // 仅处理度为2的流形边界
            return; 

        for (TopTools_ListIteratorOfListOfShape it(edgeList); it.More(); it.Next()) {
            TopoDS_Edge adjEdge = TopoDS::Edge(it.Value());
            if (!adjEdge.IsSame(shortEdge)) {
                // 如果邻边也是短边，跳过
                if (isShortEdge(adjEdge)) 
                    continue;

                double angle = GeometryUtils::ComputeIntersectionAngle(shortEdge, adjEdge, v);
                double angleDeg = angle * 180.0 / M_PI;
				//优先选择夹角大于阈值且最大的边进行合并
                if (angleDeg >= angleThresholdDeg_ && angleDeg > maxAngle) {
                    maxAngle = angleDeg;
                    bestCommonVertex = v;
                    bestAdjacentEdge = adjEdge;
                    found = true;
                }
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

    // 收集所有边
    std::vector<TopoDS_Edge> allEdges = collectAllEdges();

    // 检测短边
    for (const auto& edge : allEdges) {
        if (isShortEdge(edge)) {
            shortEdges_.push_back(edge);
        }
    }
}

TopoDS_Shape MergeShortEdge::executeUnifySameDomain(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) 
        return shape;

    // 参数：(Shape, UnifyEdges, UnifyFaces, ConcatBSplines)
    // bool参数强制合并连续的边
    ShapeUpgrade_UnifySameDomain unifier(shape, Standard_True, Standard_True, Standard_False);
    unifier.Build();

    return unifier.Shape();
}

std::vector<MergeCandidate> MergeShortEdge::findMergeCandidates()
{
    mergeCandidates_.clear();
    if (inputShapes_.empty()) 
        return mergeCandidates_;

    buildTopologyMap();
    if (shortEdges_.empty()) 
        detectShortEdges();

    TopTools_MapOfShape processedEdges;

    for (const auto& shortEdge : shortEdges_) {
        if (processedEdges.Contains(shortEdge)) 
            continue;

        TopoDS_Vertex commonVertex;
        TopoDS_Edge adjacentEdge;

        if (findBestAdjacentEdge(shortEdge, commonVertex, adjacentEdge)) {
            // 确保不重复使用邻边
            if (processedEdges.Contains(adjacentEdge)) 
                continue;

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
    }

    return mergeCandidates_;
}

int MergeShortEdge::performMerge()
{
    if (inputShapes_.empty())
        return 0;

    BRep_Builder builder;
    TopoDS_Compound initialCompound;
    builder.MakeCompound(initialCompound);
    for (const auto& shape : inputShapes_) {
        builder.Add(initialCompound, shape);
    }

    resultShape_ = executeUnifySameDomain(initialCompound);
    inputShapes_.clear();
    inputShapes_.push_back(resultShape_);
    if (bboxRatio_ > 0.0) {
        setLengthThresholdByBoundingBox(bboxRatio_);
    }

    //清理旧缓存
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();

    // 获取合并候选项
    findMergeCandidates();

    if (mergeCandidates_.empty()) {
        return 0;
    }

    // 执行合并
    BRepTools_ReShape reshaper;
    int mergeCount = 0;
    TopTools_MapOfShape removedEdges;

    for (const auto& candidate : mergeCandidates_) {
        if (removedEdges.Contains(candidate.edge1) || removedEdges.Contains(candidate.edge2)) {
            continue;
        }

        TopoDS_Edge newEdge = createMergedEdge(candidate);
        if (newEdge.IsNull()) 
            continue;

        // 设置合适的公差，由于底层几何改变为直线，必须放大公差让 ShapeFix 修正 PCurve
        builder.UpdateEdge(newEdge, candidate.length1 * 1.5);
         
        reshaper.Replace(candidate.edge2, newEdge);
        reshaper.Remove(candidate.edge1);

        removedEdges.Add(candidate.edge1);
        removedEdges.Add(candidate.edge2);
        mergeCount++;
    }

    resultShape_ = reshaper.Apply(resultShape_);

    // 严谨的拓扑修复阶段
    ShapeFix_Shape fixer(resultShape_);
    fixer.SetPrecision(lengthThreshold_ * 2.0); // 提供基于短边阈值的修复精度
    fixer.SetMaxTolerance(lengthThreshold_ * 10.0);
    fixer.FixWireTool()->FixEdgeCurvesMode() = 1; // 强制修复边曲线
    fixer.FixFaceTool()->FixMissingSeamMode() = 1;
    fixer.Perform();

    resultShape_ = fixer.Shape();

    if (!resultShape_.IsNull()) {
        inputShapes_.clear();
        inputShapes_.push_back(resultShape_);
        if (bboxRatio_ > 0.0) 
            setLengthThresholdByBoundingBox(bboxRatio_);
        shortEdges_.clear();
        mergeCandidates_.clear();
        vertexEdgeMap_.Clear();
    }

    return mergeCount;
}

TopoDS_Edge MergeShortEdge::createMergedEdge(const MergeCandidate& candidate) const
{
    TopoDS_Vertex v1 = getOtherVertex(candidate.edge1, candidate.commonVertex);
    TopoDS_Vertex v2 = getOtherVertex(candidate.edge2, candidate.commonVertex);

    if (v1.IsNull() || v2.IsNull())
        return TopoDS_Edge();

    gp_Pnt p1 = BRep_Tool::Pnt(v1);
    gp_Pnt p2 = BRep_Tool::Pnt(v2);
    gp_Pnt p_common = BRep_Tool::Pnt(candidate.commonVertex);

    TopLoc_Location loc;
    Standard_Real first, last;
    Handle(Geom_Curve) baseCurve = BRep_Tool::Curve(candidate.edge2, loc, first, last);

    // 应用位置变换到曲线
    if (!baseCurve.IsNull() && !loc.IsIdentity()) {
        baseCurve = Handle(Geom_Curve)::DownCast(baseCurve->Transformed(loc.Transformation()));
    }

    if (!baseCurve.IsNull()) {

        // 循环剥离所有的裁剪外衣，获取真正的底层数学基曲线
        Handle(Geom_Curve) actualCurve = baseCurve;
        while (actualCurve->IsKind(STANDARD_TYPE(Geom_TrimmedCurve))) {
            actualCurve = Handle(Geom_TrimmedCurve)::DownCast(actualCurve)->BasisCurve();
        }

        try {
            // 处理直线 (Geom_Line)
            if (actualCurve->IsKind(STANDARD_TYPE(Geom_Line))) {
                // 对于直线，直接传入两端的顶点生成。
                // 传入 TopoDS_Vertex 而不是 gp_Pnt 是为了完美保持拓扑连通性
                BRepBuilderAPI_MakeEdge edgeMaker(v1, v2);
                if (edgeMaker.IsDone()) {
                    TopoDS_Edge newEdge = edgeMaker.Edge();
                    newEdge.Orientation(candidate.edge2.Orientation()); // 保持长边方向
                    return newEdge;
                }
            }
            // 处理圆 (Geom_Circle)
            else if (actualCurve->IsKind(STANDARD_TYPE(Geom_Circle))) {
                // 将短边的远端点 p1 投影到主边所在的基础圆上，确保三点共圆
                GeomAPI_ProjectPointOnCurve projector(p1, actualCurve);
                gp_Pnt p1_proj = p1;
                if (projector.NbPoints() > 0) {
                    p1_proj = projector.NearestPoint();
                }

                // 使用三点构建圆弧（投影起点，原始公共连接点，长边终点）
                GC_MakeArcOfCircle arcMaker(p1_proj, p_common, p2);
                if (arcMaker.IsDone()) {
                    // 同样，传入已有的 v1, v2 顶点，绑定拓扑关系
                    BRepBuilderAPI_MakeEdge edgeMaker(arcMaker.Value(), v1, v2);
                    if (edgeMaker.IsDone()) {
                        TopoDS_Edge newEdge = edgeMaker.Edge();
                        newEdge.Orientation(candidate.edge2.Orientation());
                        return newEdge;
                    }
                }
            }
            //处理椭圆、抛物线、双曲线
            else if (actualCurve->IsKind(STANDARD_TYPE(Geom_Ellipse)) ||
                actualCurve->IsKind(STANDARD_TYPE(Geom_Parabola)) ||
                actualCurve->IsKind(STANDARD_TYPE(Geom_Hyperbola))) {

                // 复用原始数学曲线，直接传入两端顶点生成新的边
                BRepBuilderAPI_MakeEdge edgeMaker(actualCurve, v1, v2);
                if (edgeMaker.IsDone()) {
                    TopoDS_Edge newEdge = edgeMaker.Edge();
                    newEdge.Orientation(candidate.edge2.Orientation());
                    return newEdge;
                }
            }
            // TODO 处理 B 样条和贝塞尔曲线 (自由曲线)
            else if (actualCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)) ||
                actualCurve->IsKind(STANDARD_TYPE(Geom_BezierCurve))) {
                return TopoDS_Edge(); // 目前不处理复杂自由曲线的合并               
            }

        }
        catch (Standard_Failure&) {
        }
    }

    // 为了保证整个模型合并流程不中断，使用两点之间的直线拟合。
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
    return TopoDS_Edge(); // 完全失败时返回空边
}

TopoDS_Vertex MergeShortEdge::getOtherVertex(const TopoDS_Edge& edge, const TopoDS_Vertex& commonVertex) const
{
    if (edge.IsNull() || commonVertex.IsNull()) 
        return TopoDS_Vertex();
    TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
    TopoDS_Vertex v2 = TopExp::LastVertex(edge);
    if (v1.IsSame(commonVertex)) 
        return v2;
    if (v2.IsSame(commonVertex)) 
        return v1;
    return TopoDS_Vertex();
}


bool MergeShortEdge::setLengthThresholdByBoundingBox(double ratio)
{
    if (inputShapes_.empty()) {
        return false;
    }

    // 检查 ratio 有效性（必须为正数且合理）
    if (ratio <= 0.0) {
        return false;
    }

    // 限制 ratio 上限，防止阈值过大（超过1%对角线长度不合理）
    if (ratio > 0.01) {
        ratio = 0.01;
    }

    // 计算包围盒对角线长度
    double diagonal = GeometryUtils::ComputeBoundingBoxDiagonal(inputShapes_);
    
    if (diagonal <= 0.0) {
        return false;
    }

    // 设置阈值为对角线长度的比例
    lengthThreshold_ = diagonal * ratio;
    
    return true;
}
