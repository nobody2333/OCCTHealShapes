#include "MergeShortEdge.h"

#include <algorithm>
#include <cmath>

#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MergeShortEdge::MergeShortEdge(double lengthThreshold, double angleThresholdDeg)
    : lengthThreshold_(lengthThreshold)
    , angleThresholdDeg_(angleThresholdDeg)
    , angleThresholdRad_(angleThresholdDeg * M_PI / 180.0)
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
}

void MergeShortEdge::setInputShapes(const std::vector<TopoDS_Shape>& shapes)
{
    inputShapes_ = shapes;
    // 过滤掉空形状
    inputShapes_.erase(
        std::remove_if(inputShapes_.begin(), inputShapes_.end(),
            [](const TopoDS_Shape& s) { return s.IsNull(); }),
        inputShapes_.end()
    );
    // 重置结果
    shortEdges_.clear();
    mergeCandidates_.clear();
    vertexEdgeMap_.Clear();
}

void MergeShortEdge::buildTopologyMap()
{
    vertexEdgeMap_.Clear();
    
    // 收集所有形状中的边并构建映射
    for (const auto& shape : inputShapes_) {
        TopologyUtils::BuildVertexEdgeMap(shape, vertexEdgeMap_);
    }
}

std::vector<TopoDS_Edge> MergeShortEdge::collectAllEdges() const
{
    std::vector<TopoDS_Edge> edges;
    
    for (const auto& shape : inputShapes_) {
        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            edges.push_back(TopoDS::Edge(exp.Current()));
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
    return (length < lengthThreshold_ + Precision::Confusion());
}

bool MergeShortEdge::findAdjacentEdge(const TopoDS_Edge& edge, TopoDS_Vertex& commonVertex, TopoDS_Edge& adjacentEdge) const
{
    if (edge.IsNull()) {
        return false;
    }

    // 获取边的两个端点
    TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
    TopoDS_Vertex v2 = TopExp::LastVertex(edge);

    if (v1.IsNull() || v2.IsNull()) {
        return false;
    }

    // 优先选择度为2的顶点（即只属于两条边的顶点）
    auto checkVertex = [&](const TopoDS_Vertex& v) -> bool {
        if (!vertexEdgeMap_.Contains(v)) {
            return false;
        }

        const TopTools_ListOfShape& edgeList = vertexEdgeMap_.FindFromKey(v);
        
        // 检查顶点度数是否为2（即只属于两条边）
        if (edgeList.Extent() != 2) {
            return false;
        }

        // 找到另一条边
        for (TopTools_ListIteratorOfListOfShape it(edgeList); it.More(); it.Next()) {
            TopoDS_Edge e = TopoDS::Edge(it.Value());
            if (!e.IsSame(edge)) {
                commonVertex = v;
                adjacentEdge = e;
                return true;
            }
        }
        return false;
    };

    // 先检查第一个顶点
    if (checkVertex(v1)) {
        return true;
    }

    // 再检查第二个顶点
    if (checkVertex(v2)) {
        return true;
    }

    return false;
}

bool MergeShortEdge::isCollinear(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex& commonVertex) const
{
    if (e1.IsNull() || e2.IsNull() || commonVertex.IsNull()) {
        return false;
    }

    // 使用 GeometryUtils 计算夹角
    double angle = GeometryUtils::ComputeIntersectionAngle(e1, e2, commonVertex);

    if (angle < 0) {
        // 计算失败
        return false;
    }

    // 将角度转换为度数进行比较
    double angleDeg = angle * 180.0 / M_PI;

    // 夹角接近180度表示共线（允许一定的容差）
    // 例如：阈值175度，则角度需要在 [175, 180] 范围内
    return (angleDeg >= angleThresholdDeg_);
}

std::vector<TopoDS_Edge> MergeShortEdge::detectShortEdges()
{
    shortEdges_.clear();

    if (inputShapes_.empty()) {
        return shortEdges_;
    }

    // 收集所有边
    std::vector<TopoDS_Edge> allEdges = collectAllEdges();

    // 检测短边
    for (const auto& edge : allEdges) {
        if (isShortEdge(edge)) {
            shortEdges_.push_back(edge);
        }
    }

    return shortEdges_;
}

std::vector<MergeCandidate> MergeShortEdge::findMergeCandidates()
{
    mergeCandidates_.clear();

    if (inputShapes_.empty()) {
        return mergeCandidates_;
    }

    // 构建拓扑映射
    buildTopologyMap();

    // 如果没有短边，先检测一次
    if (shortEdges_.empty()) {
        detectShortEdges();
    }

    // 遍历所有短边，寻找可合并的候选
    for (const auto& shortEdge : shortEdges_) {
        TopoDS_Vertex commonVertex;
        TopoDS_Edge adjacentEdge;

        // 查找相邻边
        if (!findAdjacentEdge(shortEdge, commonVertex, adjacentEdge)) {
            continue;
        }

        // 检查是否共线
        if (!isCollinear(shortEdge, adjacentEdge, commonVertex)) {
            continue;
        }

        // 创建合并候选
        MergeCandidate candidate;
        candidate.edge1 = shortEdge;
        candidate.edge2 = adjacentEdge;
        candidate.commonVertex = commonVertex;
        candidate.length1 = GeometryUtils::ComputeEdgeLength(shortEdge);
        candidate.length2 = GeometryUtils::ComputeEdgeLength(adjacentEdge);
        candidate.angle = GeometryUtils::ComputeIntersectionAngle(shortEdge, adjacentEdge, commonVertex);

        mergeCandidates_.push_back(candidate);
    }

    return mergeCandidates_;
}

bool MergeShortEdge::performMerge()
{
    // TODO: 实现实际的边合并逻辑
    // 这需要使用 OpenCASCADE 的形状修复或建模算法
    // 例如：BRepBuilderAPI_MakeEdge, BRepAlgoAPI_Fuse 等
    
    // 目前仅返回候选数量 > 0
    if (mergeCandidates_.empty()) {
        findMergeCandidates();
    }
    
    return !mergeCandidates_.empty();
}

bool MergeShortEdge::exportStep(const std::string& fileName)
{
    // TODO: 实现 STEP 导出功能
    // 可以使用 STEPCAFControl_Writer 或 STEPControl_Writer
    (void)fileName; // 避免未使用参数警告
    return false;
}
