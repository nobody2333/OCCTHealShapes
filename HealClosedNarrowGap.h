#pragma once

#include <string>
#include <vector>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>

// 处理由 face 边界组成的闭合多边细缝。
// 第一版只处理由单面边界 edge 组成的闭合 loop，并尝试用一侧 edge chain 替换另一侧 edge chain。
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

    ChainChoice chooseChainToKeep(const ClosedGapCandidate& candidate) const;
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
    int findSameEdgeIndex(const std::vector<TopoDS_Edge>& edges, const TopoDS_Edge& edge) const;
    bool ensurePCurveOnFace(TopoDS_Edge& edge, const TopoDS_Face& face, double tolerance, std::string& errorMessage) const;
    void appendMessage(const std::string& message);
};
