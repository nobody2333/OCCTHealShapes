#pragma once

#include <string>
#include <vector>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// 修复前的保守拓扑整理：近重合顶点合并、窄边界短端直接修复。
class ShapePreProcess {
public:
    struct Params {
        bool enableDuplicateVertexMerge = true;
        bool enableShortBoundaryEdgeCollapse = false;
        // 仅供近重合顶点合并使用，不作为短边塌缩预算。
        double zeroDistanceTolerance = 1.0e-4;
        // Zero/One Cusp 端部连接边必须不超过此长度，且明显短于两侧长边。
        double maxShortBoundaryEdgeLength = 1.5e-1;
        double minLongToShortLengthRatio = 3.0;
        // 单侧旧边链间距上限；双侧新边到两条旧边界的距离上限，不是 edge 容差。
        double maxBoundaryChange = 1.5e-1;
        double maxEdgeTolerance = 1.0e-2;
        // 双侧共享投影边到另一原支撑面可允许的实测偏差；与单侧 edge 预算分开。
        double maxBilateralSharedEdgeTolerance = 1.0e-1;
        // 为 true 时，如果两个 vertex 分属不同 solid，则跳过。
        bool forbidAcrossDifferentSolids = false;
        // 为 true 时，单次替换后运行 ShapeFix_Shape；默认关闭以减少额外改动。
        bool runShapeFixAfterMerge = false;
        // 为 true 时增加少量候选跳过原因；关键提交结果始终记录。
        bool verbose = false;
    };

    explicit ShapePreProcess(const Params& params = Params());

    // 设置输入 shape，并清空上一次执行状态。
    void setInputShape(const TopoDS_Shape& shape);

    // 获取最近一次执行后的结果；未修改时等于输入 shape。
    TopoDS_Shape getResultShape() const { return resultShape_; }

    // 获取最近一次输入 shape。
    TopoDS_Shape getInputShape() const { return inputShape_; }

    // 最近一次执行是否发生修改。
    bool isModified() const { return modified_; }

    // 最近一次执行成功合并的 vertex pair 数量。
    int getMergedCount() const { return mergedCount_; }
    int getCollapsedBoundaryCount() const { return collapsedBoundaryCount_; }

    // 最近一次执行的过程日志。
    const std::vector<std::string>& getMessages() const { return messages_; }

    // 使用当前 inputShape_ 执行预处理，返回两类已提交操作数之和。
    int perform();

    // 便捷入口：先 setInputShape(shape)，再执行 perform()。
    int perform(const TopoDS_Shape& shape);

private:
    // 单个 pair 候选，只表示两个 vertex 的一次合并尝试。
    struct Candidate {
        // 保留下来的 vertex。
        TopoDS_Vertex keptVertex;
        // 被替换掉的重复 vertex。
        TopoDS_Vertex duplicateVertex;
        // 两个 vertex 的 3D 点距离。
        double distance = 0.0;
        // 候选置信度，越大越优先尝试。
        double confidence = 0.0;
    };

    struct BoundaryPath {
        std::vector<TopoDS_Edge> edges;
        std::vector<TopoDS_Vertex> vertices;
        double length = 0.0;
    };

    struct BoundaryGapCandidate {
        BoundaryPath keptPath;
        BoundaryPath movedPath;
        std::vector<TopoDS_Edge> shortEdges;
        TopoDS_Face keptFace;
        TopoDS_Face defectiveFace;
        double maxWidth = 0.0;
    };

    // 从当前 shape 中提取孤立的 near-zero vertex pair 候选。
    std::vector<Candidate> findDuplicateVertexCandidates(const TopoDS_Shape& shape);

    // 判断两个 vertex 是否满足 near-zero 距离条件。
    bool isZeroDistanceVertexPair(
        const TopoDS_Vertex& firstVertex,
        const TopoDS_Vertex& secondVertex,
        double& distance) const;

    // 根据局部拓扑和容差策略选择保留哪个 vertex。
    TopoDS_Vertex chooseKeptVertex(
        const TopoDS_Vertex& firstVertex,
        const TopoDS_Vertex& secondVertex,
        const TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap) const;

    // 检查合并 pair 是否会造成明显拓扑破坏。
    bool isSafeToMergePair(
        const TopoDS_Shape& shape,
        const Candidate& candidate,
        const TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
        std::string& errorMessage) const;

    // 尝试提交一个 pair 的合并；失败时不修改 currentShape。
    bool tryMergeCandidate(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& mergedShape,
        std::string& errorMessage);

    std::vector<BoundaryGapCandidate> findBoundaryGapCandidates(
        const TopoDS_Shape& shape) const;

    bool tryRepairBoundaryGap(
        const TopoDS_Shape& currentShape,
        const BoundaryGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage) const;

    // 对替换后的 shape 做尽量局部的后处理。
    bool repairAffectedTopology(
        TopoDS_Shape& shape,
        const Candidate& candidate,
        std::string& errorMessage) const;

    // 校验替换结果没有引入新的整体拓扑问题。
    bool validateMergeResult(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& afterShape,
        double collapsedEdgeTolerance,
        std::string& errorMessage) const;

    // 计算当前 shape 的自动迭代上限。
    int computeAutomaticIterationLimit(const TopoDS_Shape& shape) const;

    // 统计当前 shape 中首尾重合或几何长度近似为 0 的 edge 数量。
    int countCollapsedEdges(const TopoDS_Shape& shape, double tolerance) const;

    // 追加过程日志。
    void appendMessage(const std::string& message);

private:
    // 当前修复参数。
    Params params_;

    // 当前输入 shape。
    TopoDS_Shape inputShape_;

    // 最近一次执行后的结果 shape。
    TopoDS_Shape resultShape_;

    // 最近一次执行是否发生修改。
    bool modified_ = false;

    // 最近一次执行成功合并的 pair 数量。
    int mergedCount_ = 0;
    int collapsedBoundaryCount_ = 0;

    // 最近一次执行的过程日志。
    std::vector<std::string> messages_;
};
