#pragma once

#include <string>
#include <vector>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// 只处理“两个不同 vertex 的几何位置近似重合”的保守修复器。
class MergeDuplicateVertex {
public:
    // 修复参数：第一版只做 pair-only / near-zero-only，不处理三个及以上重复点聚类。
    struct Params {
        // 两点距离小于该值时才认为是近似 0 距离。
        double zeroDistanceTolerance = 1.0e-7;
        // 为 true 时，如果两个 vertex 分属不同 solid，则跳过。
        bool forbidAcrossDifferentSolids = true;
        // 为 true 时，单次替换后运行 ShapeFix_Shape；默认关闭以减少额外改动。
        bool runShapeFixAfterMerge = false;
        // 预留调试开关；当前 messages_ 始终保留关键过程日志。
        bool verbose = false;
    };

    explicit MergeDuplicateVertex(const Params& params = Params());

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

    // 最近一次执行的过程日志。
    const std::vector<std::string>& getMessages() const { return messages_; }

    // 使用当前 inputShape_ 执行修复，返回成功合并的 pair 数量。
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

    // 最近一次执行的过程日志。
    std::vector<std::string> messages_;
};
