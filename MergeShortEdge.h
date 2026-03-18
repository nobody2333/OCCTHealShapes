#pragma once

#include <vector>
#include <memory>

#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// 合并候选结构
struct MergeCandidate {
    TopoDS_Edge edge1;           // 第一条边（通常是短边）
    TopoDS_Edge edge2;           // 第二条边（与 edge1 相邻）
    TopoDS_Vertex commonVertex;  // 公共顶点
    double length1;              // edge1 的长度
    double length2;              // edge2 的长度
    double angle;                // 两边夹角（弧度）

    MergeCandidate() : length1(0.0), length2(0.0), angle(0.0) {}
};

/**
 * @brief 短边检测与合并类
 * 
 * 使用 GeometryUtils 计算边长和夹角
 * 使用 TopologyUtils 构建顶点到边的拓扑映射
 */
class MergeShortEdge {
public:
    /**
     * @brief 默认构造函数
     * @param lengthThreshold 短边长度阈值（默认 0.5）
     * @param angleThreshold  共线角度阈值（默认 175 度，转换为弧度）
     */
    MergeShortEdge(double lengthThreshold = 0.5, double angleThresholdDeg = 175.0);
    
    ~MergeShortEdge() = default;

    // 禁止拷贝，允许移动
    MergeShortEdge(const MergeShortEdge&) = delete;
    MergeShortEdge& operator=(const MergeShortEdge&) = delete;
    MergeShortEdge(MergeShortEdge&&) = default;
    MergeShortEdge& operator=(MergeShortEdge&&) = default;

    /**
     * @brief 设置输入形状（可以是 Solid、Shell、Face 或 Compound）
     * @param shape 输入的 TopoDS_Shape
     */
    void setInputShape(const TopoDS_Shape& shape);

    /**
     * @brief 设置多个输入形状
     * @param shapes 输入形状列表
     */
    void setInputShapes(const std::vector<TopoDS_Shape>& shapes);

    // 设置阈值参数
    void setLengthThreshold(double threshold) { lengthThreshold_ = threshold; }
    void setAngleThresholdDeg(double degrees) { angleThresholdDeg_ = degrees; }
    
    // 获取当前阈值
    double getLengthThreshold() const { return lengthThreshold_; }
    double getAngleThresholdDeg() const { return angleThresholdDeg_; }

    /**
     * @brief 检测所有短边
     * @return 检测到的短边列表
     */
    std::vector<TopoDS_Edge> detectShortEdges();

    /**
     * @brief 获取可合并的短边候选对
     * @return 满足条件的合并候选列表
     */
    std::vector<MergeCandidate> findMergeCandidates();

    /**
     * @brief 执行短边合并（待实现）
     * @return 是否成功
     */
    bool performMerge();

    // 获取结果
    const std::vector<TopoDS_Edge>& getShortEdges() const { return shortEdges_; }
    const std::vector<MergeCandidate>& getMergeCandidates() const { return mergeCandidates_; }
    const TopTools_IndexedDataMapOfShapeListOfShape& getVertexEdgeMap() const { return vertexEdgeMap_; }

    /**
     * @brief 导出 STEP 文件（待实现）
     * @param fileName 输出文件名
     * @return 是否成功
     */
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

    // 内部辅助方法
    void buildTopologyMap();  // 构建顶点到边的映射
    
    /**
     * @brief 检查边是否为短边
     * @param edge 待检查的边
     * @return 如果是短边返回 true
     */
    bool isShortEdge(const TopoDS_Edge& edge) const;

    /**
     * @brief 检查边在公共顶点处的连接性
     * @param edge 待检查的边
     * @param commonVertex 输出公共顶点
     * @param adjacentEdge 输出相邻边
     * @return 如果找到相邻边返回 true
     */
    bool findAdjacentEdge(const TopoDS_Edge& edge, TopoDS_Vertex& commonVertex, TopoDS_Edge& adjacentEdge) const;

    /**
     * @brief 检查两边是否满足共线条件（夹角接近180度）
     * @param e1 第一条边
     * @param e2 第二条边
     * @param commonVertex 公共顶点
     * @return 如果共线返回 true
     */
    bool isCollinear(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex& commonVertex) const;

    /**
     * @brief 收集形状中的所有边
     * @return 边的列表
     */
    std::vector<TopoDS_Edge> collectAllEdges() const;
};
