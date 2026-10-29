#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "RepairDiagnostics.h"

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>

class BRepTools_ReShape;

// 处理由 face 边界组成的闭合多边细缝。
// 支持尖点端部和短 edge 连接端部，并优先尝试用中线 edge chain 收缩两侧边界。
class HealClosedNarrowGap {
public:
    // 指定中线 edge chain 的构造方式；Auto 保持原有自动分派行为。
    enum class MidlineConstructionMode {
        Auto,
        LegacyOneToOne,
        ProjectionCorrespondence
    };

    struct Params {
        // 两侧 chain 之间允许的最大宽度。
        double maxGapWidth = 0.8;
        // 候选两侧 chain 长度相似度下限：min(lengthA, lengthB) / max(lengthA, lengthB)。
        double minLengthRatio = 0.6;
        // 候选细缝的最小长宽比；值越大越偏向狭长细缝。
        double minLengthWidthRatio = 5.0;
        // 识别细缝端部尖点时允许的最大夹角。
        double maxTipAngle = 35.0 * 3.14159265358979323846 / 180.0;
        // 宽度检测的目标采样间距，单位与模型单位一致。
        // 内部自适应取 20～120 点；间距或长度无效时回退到 20 点。
        double targetWidthSampleSpacing = 0.5;
        // 中线曲线的目标采样间距，单位与模型单位一致。
        // 内部自适应取 8～80 点；间距或长度无效时回退到 20 点。
        double targetSampleSpacing = 1.0;
        // 中线构造模式；投影对应模式也可以显式用于拓扑分段对称的候选。
        MidlineConstructionMode midlineConstructionMode =
            MidlineConstructionMode::Auto;
        // 投影点靠近对侧已有顶点时，允许吸附成真实顶点配对的几何容差。
        double projectionVertexSnapTolerance = 1.0e-4;
        // 中线收缩失败时，是否允许删除 wire 异常的 face，并用剩余边界拟合新 face。
        bool useDeleteAndFillDefectiveFace = true;
        // DeleteAndFillDefectiveFace 拟合 B-spline 支撑面的工作容差。
        double deleteAndFillTolerance = 1.0e-3;
        // DeleteAndFillDefectiveFace 允许的最大边界 G0 偏差。
        double deleteAndFillMaxG0Error = 1.0e-1;
        // 中线构造、face 重建及最终 ShapeFix_Shape 使用的修复容差。
        double midlineRepairTolerance = 0.1;
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
    // 闭合细缝端部可以是两侧共用尖点，也可以是连接两侧的短 edge。
    enum class GapEndType {
        SharedTip,
        ShortEdge
    };

    struct GapEnd {
        GapEndType type = GapEndType::SharedTip;
        TopoDS_Vertex vertexA;
        TopoDS_Vertex vertexB;
        TopoDS_Edge connectorEdge;
        TopoDS_Face connectorOwnerFace;
        double score = 0.0;
    };

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
        GapEnd startEnd;
        GapEnd finishEnd;
        double meanWidth = 0.0;
        double maxWidth = 0.0;
        double maxWidthAToB = 0.0;
        double maxWidthBToA = 0.0;
        int usedWidthSampleCount = 0;
        double lengthRatio = 0.0;
        double lengthWidthRatio = 0.0;
        double endpointSpan = 0.0;
        double endpointSpanWidthRatio = 0.0;
        double endpointScore = 0.0;
        double confidence = 0.0;
    };

    struct ChainChoice {
        EdgeChain keptChain;
        EdgeChain replaceChain;
    };

    // 保存一侧 chain 按真实 edge 弧长计算的归一化顶点站点。
    struct ChainStationData {
        std::vector<double> edgeLengths;
        std::vector<double> vertexStations;
        double totalLength = 0.0;
    };

    // 一个源点投影到整条目标 chain 后的最近解及吸附信息。
    struct ChainProjectionResult {
        gp_Pnt projectedPoint;
        gp_Pnt rawProjectedPoint;
        std::size_t targetEdgeIndex =
            std::numeric_limits<std::size_t>::max();
        double targetParameter = 0.0;
        double targetStation = 0.0;
        double rawTargetStation = 0.0;
        double distance = 0.0;
        double rawDistance = 0.0;
        std::size_t snappedTargetVertexIndex =
            std::numeric_limits<std::size_t>::max();
        bool snappedToVertex = false;
    };

    // 表示 chainA(stationA) 与 chainB(stationB) 的一个单调对应节点。
    struct ProjectionCorrespondenceNode {
        double stationA = 0.0;
        double stationB = 0.0;
        gp_Pnt pointA;
        gp_Pnt pointB;
        std::size_t sourceVertexIndexA =
            std::numeric_limits<std::size_t>::max();
        std::size_t sourceVertexIndexB =
            std::numeric_limits<std::size_t>::max();
        double projectionDistance = 0.0;
        bool hasSourceVertexA = false;
        bool hasSourceVertexB = false;
        bool snappedPair = false;
        // 该节点是否由目标短边检查合并得到；它与投影吸附配对是两种不同机制。
        bool mergedForShortEdge = false;
        bool endpoint = false;
    };

    // 双向投影中经过确定性冲突消解后保留的一对真实顶点。
    struct SnappedVertexPair {
        std::size_t sourceVertexIndexA =
            std::numeric_limits<std::size_t>::max();
        std::size_t sourceVertexIndexB =
            std::numeric_limits<std::size_t>::max();
        double distance = 0.0;
    };

    // 中线构造结果；映射中保存规范中线 edge 的索引，避免各 face 持有不同副本。
    struct MidlineChainData {
        EdgeChain chain;
        std::vector<std::vector<std::size_t>> replacementsForChainA;
        std::vector<std::vector<std::size_t>> replacementsForChainB;
        // A、B 两侧源顶点数量可能不同，必须分别映射到统一目标顶点。
        std::vector<std::size_t> targetVertexIndicesForChainA;
        std::vector<std::size_t> targetVertexIndicesForChainB;
        std::vector<ProjectionCorrespondenceNode> correspondenceNodes;
        int shortEdgeMergeCount = 0;
        MidlineConstructionMode constructionMode =
            MidlineConstructionMode::LegacyOneToOne;
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
        // 非尖点端部收缩后，短连接 edge 必须从所属 wire 中移除。
        std::vector<TopoDS_Edge> edgesToRemove;
    };

    // 删除缺陷 face 后，用仍被其他 face 使用的边和细缝另一侧边组成新的补面边界。
    struct DeleteAndFillPlan {
        TopoDS_Face defectiveFace;
        EdgeChain boundary;
        std::vector<TopoDS_Edge> removedEdges;
        int missingWireEdgeUses = 0;
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
        const RepairDiagnostics::FailureList& existingFailures,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool tryRepairByMidlineCollapse(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool tryRepairByDeletingAndFillingDefectiveFace(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool buildDeleteAndFillPlan(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        DeleteAndFillPlan& plan,
        std::string& errorMessage) const;

    bool buildDeleteAndFillPatchFace(
        const DeleteAndFillPlan& plan,
        TopoDS_Face& patchFace,
        double& g0Error,
        std::string& errorMessage) const;

    // 根据保留边界另一侧 face 的 edge-use 方向确定补面方向，避免继承异常旧 face 的法线。
    bool orientDeleteAndFillPatchFace(
        const DeleteAndFillPlan& plan,
        TopoDS_Face& patchFace,
        bool& reversedPatchFace,
        std::string& errorMessage) const;

    // 补面提交前检查每条边界 edge 在两个相邻 face 中是否以相反方向使用。
    bool validateDeleteAndFillBoundaryOrientations(
        const TopoDS_Shape& shape,
        const EdgeChain& boundary,
        std::string& errorMessage) const;

    bool tryRepairBySideReplacement(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const ClosedGapCandidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    ChainChoice chooseChainToKeep(const ClosedGapCandidate& candidate) const;
    bool sampleChainPoint(const EdgeChain& chain, double normalized, gp_Pnt& point) const;

    // 在整条实际裁剪 chain 上查找离 sourcePoint 最近的点。
    bool findNearestPointOnChain(
        const gp_Pnt& sourcePoint,
        const EdgeChain& targetChain,
        std::size_t& targetEdgeIndex,
        double& targetParameter,
        gp_Pnt& targetPoint,
        double& distance,
        std::string& errorMessage) const;

    // 两侧分别采样并投影到对侧 chain，计算离散双向最近距离。
    bool evaluateBidirectionalGap(
        const EdgeChain& chainA,
        const EdgeChain& chainB,
        double& maxWidthAToB,
        double& maxWidthBToA,
        double& meanWidth,
        int& usedSampleCount,
        std::string& errorMessage) const;

    bool buildMidlineChain(
        const ClosedGapCandidate& candidate,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    bool isCandidateEndConsistent(
        const GapEnd& end,
        const TopoDS_Vertex& vertexA,
        const TopoDS_Vertex& vertexB) const;

    // 尖点端部复用原 vertex；短 edge 端部在两个端点的几何中点创建统一 vertex。
    bool buildTargetEndVertex(
        const GapEnd& end,
        TopoDS_Vertex& targetVertex,
        std::string& errorMessage) const;

    bool canUseLegacyOneToOnePath(
        const ClosedGapCandidate& candidate) const;

    bool resolveMidlineConstructionMode(
        const ClosedGapCandidate& candidate,
        MidlineConstructionMode& mode,
        std::string& errorMessage) const;

    bool buildMidlineChainOneToOne(
        const ClosedGapCandidate& candidate,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    bool buildMidlineChainByProjectionCorrespondence(
        const ClosedGapCandidate& candidate,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    bool buildNormalizedVertexStations(
        const EdgeChain& chain,
        const char* chainName,
        ChainStationData& stationData,
        std::string& errorMessage) const;

    bool projectPointToChain(
        const gp_Pnt& sourcePoint,
        const EdgeChain& targetChain,
        const ChainStationData& targetStationData,
        ChainProjectionResult& projection,
        std::string& errorMessage) const;

    bool selectSnappedVertexPairs(
        const EdgeChain& chainA,
        const EdgeChain& chainB,
        const std::vector<ChainProjectionResult>& projectionsFromA,
        const std::vector<ChainProjectionResult>& projectionsFromB,
        std::vector<SnappedVertexPair>& selectedPairs,
        std::string& errorMessage) const;

    bool buildProjectionCorrespondenceNodes(
        const ClosedGapCandidate& candidate,
        const ChainStationData& chainAStations,
        const ChainStationData& chainBStations,
        std::vector<ProjectionCorrespondenceNode>& nodes,
        std::string& errorMessage) const;

    bool mergeProjectionCorrespondenceNodes(
        std::vector<ProjectionCorrespondenceNode>& nodes,
        const EdgeChain& chainA,
        const EdgeChain& chainB,
        const ChainStationData& chainAStations,
        const ChainStationData& chainBStations,
        std::string& errorMessage) const;

    bool buildProjectionTargetVertices(
        const ClosedGapCandidate& candidate,
        const std::vector<ProjectionCorrespondenceNode>& nodes,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    bool buildProjectionTargetEdges(
        const ClosedGapCandidate& candidate,
        const ChainStationData& chainAStations,
        const ChainStationData& chainBStations,
        const std::vector<ProjectionCorrespondenceNode>& nodes,
        MidlineChainData& midlineData,
        std::string& errorMessage) const;

    // 根据已生成目标 edge 的实际长度和端点容差，合并可能在后续修复中退化的相邻节点。
    bool mergeProjectionNodesForShortTargetEdges(
        std::vector<ProjectionCorrespondenceNode>& nodes,
        const MidlineChainData& midlineData,
        bool& nodesMerged,
        std::string& errorMessage) const;

    bool buildSourceToTargetEdgeMapping(
        const std::vector<double>& sourceVertexStations,
        const std::vector<double>& targetStations,
        const char* chainName,
        std::vector<std::vector<std::size_t>>& mapping,
        std::string& errorMessage) const;

    bool validateReplacementMapping(
        const std::vector<std::vector<std::size_t>>& mapping,
        std::size_t sourceEdgeCount,
        std::size_t targetEdgeCount,
        const char* chainName,
        std::string& errorMessage) const;

    static const char* midlineConstructionModeName(
        MidlineConstructionMode mode);

    static const char* gapEndTypeName(GapEndType type);

    bool rebuildFacesByCollapsingToMidline(
        const TopoDS_Shape& currentShape,
        const ClosedGapCandidate& candidate,
        const MidlineChainData& midlineData,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    bool rebuildFaceByReplacingMappedEdges(
        const TopoDS_Face& faceToRebuild,
        const std::vector<EdgeReplacementTask>& edgeTasks,
        const std::vector<TopoDS_Edge>& edgesToRemove,
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
        const std::vector<TopoDS_Edge>& removedEdges,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    // 本轮候选复用 currentShape 的错误基线；中间 trial 比较仍使用原重载。
    bool validateRepairResult(
        const RepairDiagnostics::FailureList& existingFailures,
        const TopoDS_Shape& afterShape,
        const EdgeChain& keptChain,
        const std::vector<TopoDS_Edge>& removedEdges,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    int computeAutomaticIterationLimit(const TopoDS_Shape& shape) const;
    // 根据两侧 chain 长度自适应取 20～120 点；间距或长度无效时回退到 20 点。
    int computeWidthSampleCount(double lengthA, double lengthB) const;
    int computeMidlineSampleCount(double lengthA, double lengthB) const;
    int findSameEdgeIndex(const std::vector<TopoDS_Edge>& edges, const TopoDS_Edge& edge) const;
    bool ensurePCurveOnFace(
        TopoDS_Edge& edge,
        const TopoDS_Face& face,
        double tolerance,
        std::string& errorMessage) const;
    void appendMessage(const std::string& message);
};
