#pragma once

#include <vector>
#include <memory>

#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

/// <summary>
/// 合并候选结构，存储两条相邻边的信息
/// </summary>
struct MergeCandidate {
    TopoDS_Edge edge1;           // 第一条边（通常是短边）
    TopoDS_Edge edge2;           // 第二条边（与 edge1 相邻）
    TopoDS_Vertex commonVertex;  // 公共顶点
    double length1;              // edge1 的长度
    double length2;              // edge2 的长度
    double angle;                // 两边夹角（弧度）

    MergeCandidate() : length1(0.0), length2(0.0), angle(0.0) {}
};

/// <summary>
/// 短边检测与合并类
/// 使用 GeometryUtils 计算边长和夹角
/// 使用 TopologyUtils 构建顶点到边的拓扑映射
/// </summary>
class MergeShortEdge {
public:
    /// <summary>
    /// 默认构造函数
    /// </summary>
    /// <param name="lengthThreshold">短边长度阈值（默认 0.5）</param>
    /// <param name="angleThresholdDeg">共线角度阈值（默认 175 度）</param>
    MergeShortEdge(double lengthThreshold = 0.5, double angleThresholdDeg = 175.0);
    
    ~MergeShortEdge() = default;

    // 禁止拷贝，允许移动
    MergeShortEdge(const MergeShortEdge&) = delete;
    MergeShortEdge& operator=(const MergeShortEdge&) = delete;
    MergeShortEdge(MergeShortEdge&&) = default;
    MergeShortEdge& operator=(MergeShortEdge&&) = default;

    /// <summary>
    /// 设置输入形状（可以是 Solid、Shell、Face 或 Compound）
    /// </summary>
    /// <param name="shape">输入的 TopoDS_Shape</param>
    void setInputShape(const TopoDS_Shape& shape);

    /// <summary>
    /// 设置多个输入形状
    /// </summary>
    /// <param name="shapes">输入形状列表</param>
    void setInputShapes(const std::vector<TopoDS_Shape>& shapes);

    /// <summary>
    /// 设置长度阈值
    /// </summary>
    /// <param name="threshold">长度阈值</param>
    void setLengthThreshold(double threshold) { lengthThreshold_ = threshold; }

    /// <summary>
    /// 设置角度阈值（度数）
    /// </summary>
    /// <param name="degrees">角度阈值（度）</param>
    void setAngleThresholdDeg(double degrees) { angleThresholdDeg_ = degrees; }
    
    /// <summary>
    /// 获取当前长度阈值
    /// </summary>
    /// <returns>长度阈值</returns>
    double getLengthThreshold() const { return lengthThreshold_; }

    /// <summary>
    /// 获取当前角度阈值（度数）
    /// </summary>
    /// <returns>角度阈值（度）</returns>
    double getAngleThresholdDeg() const { return angleThresholdDeg_; }

    /// <summary>
    /// 检测所有短边
    /// </summary>
    /// <returns>检测到的短边列表</returns>
    std::vector<TopoDS_Edge> detectShortEdges();

    /// <summary>
    /// 获取可合并的短边候选对
    /// </summary>
    /// <returns>满足条件的合并候选列表</returns>
    std::vector<MergeCandidate> findMergeCandidates();

    /// <summary>
    /// 执行短边合并（TODO）
    /// </summary>
    /// <returns>是否成功</returns>
    bool performMerge();

    /// <summary>
    /// 获取检测到的短边列表
    /// </summary>
    /// <returns>短边列表</returns>
    const std::vector<TopoDS_Edge>& getShortEdges() const { return shortEdges_; }

    /// <summary>
    /// 获取可合并候选列表
    /// </summary>
    /// <returns>合并候选列表</returns>
    const std::vector<MergeCandidate>& getMergeCandidates() const { return mergeCandidates_; }

    /// <summary>
    /// 获取顶点到边的拓扑映射
    /// </summary>
    /// <returns>拓扑映射</returns>
    const TopTools_IndexedDataMapOfShapeListOfShape& getVertexEdgeMap() const { return vertexEdgeMap_; }

    /// <summary>
    /// 导出 STEP 文件（待实现）
    /// </summary>
    /// <param name="fileName">输出文件名</param>
    /// <returns>是否成功</returns>
    bool exportStep(const std::string& fileName);

private:
    // 阈值参数
    double lengthThreshold_;      // 短边长度阈值
    double angleThresholdDeg_;    // 共线角度阈值（度）
    double angleThresholdRad_;    // 共线角度阈值（弧度）

    // 输入数据
    std::vector<TopoDS_Shape> inputShapes_;  // 输入形状列表
    
    // 拓扑映射
    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap_;  // 顶点到边的映射
    
    // 结果数据
    std::vector<TopoDS_Edge> shortEdges_;       // 检测到的短边
    std::vector<MergeCandidate> mergeCandidates_;  // 可合并的候选

    /// <summary>
    /// 构建顶点到边的映射
    /// </summary>
    void buildTopologyMap();
    
    /// <summary>
    /// 检查边是否为短边
    /// </summary>
    /// <param name="edge">待检查的边</param>
    /// <returns>如果是短边返回 true</returns>
    bool isShortEdge(const TopoDS_Edge& edge) const;

    /// <summary>
    /// 检查边在公共顶点处的连接性，查找相邻边
    /// </summary>
    /// <param name="edge">待检查的边</param>
    /// <param name="commonVertex">输出公共顶点</param>
    /// <param name="adjacentEdge">输出相邻边</param>
    /// <returns>如果找到相邻边返回 true</returns>
    bool findAdjacentEdge(const TopoDS_Edge& edge, TopoDS_Vertex& commonVertex, TopoDS_Edge& adjacentEdge) const;

    /// <summary>
    /// 检查两边是否满足共线条件（夹角接近180度）
    /// </summary>
    /// <param name="e1">第一条边</param>
    /// <param name="e2">第二条边</param>
    /// <param name="commonVertex">公共顶点</param>
    /// <returns>如果共线返回 true</returns>
    bool isCollinear(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex& commonVertex) const;

    /// <summary>
    /// 收集形状中的所有边
    /// </summary>
    /// <returns>边的列表</returns>
    std::vector<TopoDS_Edge> collectAllEdges() const;
};
