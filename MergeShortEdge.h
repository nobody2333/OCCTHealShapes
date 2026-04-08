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
/// 可选参数包括包围盒对角线比例和两边夹角
/// 长度阈值可指定具体值，也可以基于包围盒对角线的比例自动计算
/// </summary>
class MergeShortEdge {
public:
    /// <summary>
    /// 默认构造函数
    /// </summary>
    /// <param name="bboxRatio">短边阈值相对于包围盒对角线的比例（默认 0.001，即千分之一）</param>
    /// <param name="angleThresholdDeg">共线角度阈值（默认 175 度）</param>
    MergeShortEdge(double bboxRatio = 0.001, double angleThresholdDeg = 175.0);
    ~MergeShortEdge() = default;

    // 禁止拷贝，允许移动
    MergeShortEdge(const MergeShortEdge&) = delete;
    MergeShortEdge& operator=(const MergeShortEdge&) = delete;
    MergeShortEdge(MergeShortEdge&&) = default;
    MergeShortEdge& operator=(MergeShortEdge&&) = default;

    /// <summary>
    /// 设置输入形状
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
    /// 基于模型包围盒对角线设置长度阈值
    /// 短边阈值 = 包围盒对角线长度 × ratio
    /// </summary>
    /// <param name="ratio">比例系数，默认0.001（千分之一）</param>
    /// <returns>是否成功设置（需要已设置输入形状）</returns>
    bool setLengthThresholdByBoundingBox(double ratio = 0.001);

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
    std::vector<TopoDS_Edge> getShortEdges();

    /// <summary>
    /// 获取可合并的短边候选对
    /// </summary>
    /// <returns>满足条件的合并候选列表</returns>
    std::vector<MergeCandidate> findMergeCandidates();

    /// <summary>
    /// 执行短边合并
    /// </summary>
    /// <returns>成功合并的边对数量</returns>
    int performMerge();

    /// <summary>
    /// 获取合并后的结果形状
    /// </summary>
    /// <returns>处理后的形状</returns>
    TopoDS_Shape getResultShape() const { return resultShape_; }

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

private:
    double lengthThreshold_;      // 短边长度阈值（实际计算值）
    double angleThresholdDeg_;    // 共线角度阈值（度）
    double bboxRatio_;            // 包围盒对角线比例（用于计算lengthThreshold_）

    std::vector<TopoDS_Shape> inputShapes_;  // 输入形状列表
    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap_;  // 顶点到边的映射
    
    std::vector<TopoDS_Edge> shortEdges_;       // 检测到的短边
    std::vector<MergeCandidate> mergeCandidates_;  // 可合并的候选
    TopoDS_Shape resultShape_;                  // 合并后的结果形状

    TopoDS_Shape executeUnifySameDomain(const TopoDS_Shape& shape);
    /// <summary>
	/// 检测所有短边，更新 shortEdges_ 列表
    /// </summary>
    /// <returns></returns>
    void detectShortEdges();

    /// <summary>
    /// 创建合并后的新边，连接两条边的非公共端点
    /// </summary>
    /// <param name="candidate">合并候选</param>
    /// <returns>新创建的边，失败返回空边</returns>
    TopoDS_Edge createMergedEdge(const MergeCandidate& candidate) const;

    /// <summary>
    /// 获取边的另一个端点（非指定顶点）
    /// </summary>
    /// <param name="edge">边</param>
    /// <param name="commonVertex">已知端点</param>
    /// <returns>另一个端点</returns>
    TopoDS_Vertex getOtherVertex(const TopoDS_Edge& edge, const TopoDS_Vertex& commonVertex) const;

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
    bool findBestAdjacentEdge(const TopoDS_Edge& shortEdge, TopoDS_Vertex& bestCommonVertex, TopoDS_Edge& bestAdjacentEdge) const;


    /// <summary>
    /// 收集形状中的所有边
    /// </summary>
    /// <returns>边的列表</returns>
    std::vector<TopoDS_Edge> collectAllEdges() const;
};
