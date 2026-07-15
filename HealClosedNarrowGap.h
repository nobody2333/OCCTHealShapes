#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>

class BRepTools_ReShape;

// 处理由 face 边界组成的闭合多边细缝。
// 第一版识别由单面边界 edge 组成的闭合 loop，并优先尝试用中线 edge chain 收缩两侧边界。
class HealClosedNarrowGap {
public:
    struct Params {
        // 两侧 chain 采样点之间允许的最大宽度。
        double maxGapWidth = 1.0e-1;
        // 候选两侧 chain 长度相似度下限：min(lengthA, lengthB) / max(lengthA, lengthB)。
        double minLengthRatio = 0.6;
        // 候选细缝的最小长宽比；值越大越偏向狭长细缝。
        double minLengthWidthRatio = 5.0;
        // 识别细缝端部尖点时允许的最大夹角。
        double maxTipAngle = 35.0 * 3.14159265358979323846 / 180.0;
        // 对两侧 chain 做同步采样的点数。
        int sampleCount = 20;
        // 是否根据每对 source edge 的长度自适应计算中线曲线采样点数。
        bool useAdaptiveSampleCount = true;
        // 中线曲线拟合使用的最小采样点数，至少包含两个端点。
        int minSampleCount = 8;
        // 中线曲线拟合使用的最大采样点数，防止长 edge 导致计算量过大。
        int maxSampleCount = 80;
        // 中线曲线的目标采样间距，单位与模型单位一致。
        double targetSampleSpacing = 1.0;
        // 为 true 时，优先用两侧采样点中线构造新 edge chain，同时替换两侧 chain。
        bool useMidlineCollapse = true;
        // 中线收缩失败时，是否回退到“一侧 chain 替换另一侧 chain”的保守方案。
        bool fallbackToSideReplacement = true;
        // 中线构造、face 重建及最终 ShapeFix_Shape 使用的修复容差。
        double midlineRepairTolerance = 1.0e-1;
        // 为 true 时优先保留总长度较长的 chain；否则优先保留 edge 数更少的 chain。
        bool preferLongerChainToKeep = true;
        // 为 true 时，替换后运行 ShapeFix_Shape 做整体后处理。
        bool runShapeFixAfterRepair = true;
        // 为 true 时，如果 free edge 数量增加，则拒绝本次修复。
        bool rejectFreeEdgeIncrease = true;
        // 为 true 时，如果 BRepCheck 出现新问题，则拒绝本次修复。
        bool rejectInvalidResult = true;
        // 调试开关；当前关键过程始终写入 messages_。
        bool verbose = false;
    };

    explicit HealClosedNarrowGap(const Params& params = Params());

    // 设置输入 shape，并清空上一轮执行状态。
    void setInputShape(const TopoDS_Shape& shape);

    // 获取最近一次执行后的结果 shape；没有成功修复时等于输入 shape。
    TopoDS_Shape getResultShape() const { return resultShape_; }

    // 获取最近一次设置的输入 shape。
    TopoDS_Shape getInputShape() const { return inputShape_; }

    // 最近一次执行是否成功修改了 shape。
    bool isModified() const { return modified_; }

    // 最近一次执行成功修复的闭合细缝数量。
    int getRepairedCount() const { return repairedCount_; }

    // 最近一次执行过程日志。
    const std::vector<std::string>& getMessages() const { return messages_; }

    // 使用当前 inputShape_ 执行修复，返回成功提交的修复次数。
    int perform();

    // 便捷入口：先 setInputShape(shape)，再 perform()。
    int perform(const TopoDS_Shape& shape);

private:
    // 细缝边界的一侧连续 edge chain，edges 和 vertices 均按采样方向排列。
    struct EdgeChain {
        std::vector<TopoDS_Edge> edges;
        std::vector<TopoDS_Vertex> vertices;
        std::vector<TopoDS_Face> ownerFaces;
        double length = 0.0;
    };

    // 一个闭合多边细缝候选。
    struct ClosedGapCandidate {
        EdgeChain chainA;
        EdgeChain chainB;
        std::vector<TopoDS_Edge> loopEdges;
        std::vector<TopoDS_Vertex> tipVertices;
        double meanWidth = 0.0;
        double maxWidth = 0.0;
        double lengthRatio = 0.0;
        double lengthWidthRatio = 0.0;
        double endpointSharpness = 0.0;
        double confidence = 0.0;
    };

    struct ChainChoice {
        EdgeChain keptChain;
        EdgeChain replaceChain;
    };

    // 中线构造结果；映射中保存规范中线 edge 的索引，避免各 face 持有不同副本。
    struct MidlineChainData {
        EdgeChain chain;
        std::vector<std::vector<std::size_t>> replacementsForChainA;
        std::vector<std::vector<std::size_t>> replacementsForChainB;
    };

    // 一个 face 内单条原 edge 的中线替换任务，替换边通过规范中线索引引用。
    struct EdgeReplacementTask {
        TopoDS_Edge sourceEdge;
        TopoDS_Vertex sourceStart;
        TopoDS_Vertex sourceEnd;
        std::vector<std::size_t> replacementEdgeIndices;
    };

    // 同一个 face 的所有局部替换必须集中后一次性重建。
    struct FaceReplacementPlan {
        TopoDS_Face face;
        std::vector<EdgeReplacementTask> edgeTasks;
    };

    Params params_;
    TopoDS_Shape inputShape_;
    TopoDS_Shape resultShape_;
    bool modified_ = false;
    int repairedCount_ = 0;
    std::vector<std::string> messages_;

    std::vector<ClosedGapCandidate> extractCandidates(const TopoDS_Shape& shape) const;
    bool tryRepairCandidate(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool tryRepairByMidlineCollapse(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool tryRepairBySideReplacement(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    ChainChoice chooseChainToKeep(const ClosedGapCandidate& candidate) const;
    bool sampleChainPoint(const EdgeChain& chain, double normalized, gp_Pnt& point) const;
    bool buildMidlineChain(
        const ClosedGapCandidate& candidate,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    bool rebuildFacesByCollapsingToMidline(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        const MidlineChainData& midlineData,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool rebuildFaceByReplacingMappedEdges(
        const TopoDS_Face& faceToRebuild,
        const std::vector<EdgeReplacementTask>& edgeTasks,
        const EdgeChain& midlineChain,
        BRepTools_ReShape& topologyReshaper,
        double tolerance,
        TopoDS_Face& rebuiltFace,
        std::string& errorMessage);

    bool rebuildFacesByReplacingChain(
        const TopoDS_Shape& currentShape,
        const EdgeChain& keptChain,
        const EdgeChain& replaceChain,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage) const;

    bool rebuildFaceByReplacingChainEdges(
        const TopoDS_Face& faceToRebuild,
        const EdgeChain& keptChain,
        const EdgeChain& replaceChain,
        double tolerance,
        TopoDS_Face& rebuiltFace,
        std::string& errorMessage) const;

    bool validateRepairResult(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& afterShape,
        const EdgeChain& keptChain,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    int computeAutomaticIterationLimit(const TopoDS_Shape& shape) const;
    // 根据一对 source edge 的长度计算中线曲线拟合点数。
    int computeMidlineSampleCount(double lengthA, double lengthB) const;
    int findSameEdgeIndex(const std::vector<TopoDS_Edge>& edges, const TopoDS_Edge& edge) const;
    bool ensurePCurveOnFace(TopoDS_Edge& edge, const TopoDS_Face& face, double tolerance, std::string& errorMessage) const;
    void appendMessage(const std::string& message);
};
