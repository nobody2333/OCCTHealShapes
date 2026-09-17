#pragma once

#include <string>
#include <vector>

#include "RepairDiagnostics.h"

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// 专门处理“同端点双边细缝”的局部细缝修复。
class HealNarrowGap {
public:
    // 修复策略：Auto 先运行常规策略；全部常规候选均失败后，再单独扫描 crossed-face 候选。
    enum class RepairStrategy {
        ReplaceOneSide,
        SurfaceIntersection,
        MidlineCollapse,
        Auto
    };

    // 修复参数：只约束本类负责的窄缝场景，不做全局 sewing 或补面。
    struct Params {
        // Auto 依次尝试单边替换、中线塌缩，再执行 crossed-face 后置筛选。
        RepairStrategy repairStrategy = RepairStrategy::Auto;
        // 是否在 Auto 的常规策略全部停止后扫描 crossed-face 候选。
        bool enableCrossedFaceBoundaryRepair = false;
        // 候选搜索阈值：用于判断两条 edge 的双向投影距离是否足够近。
        double candidateMaxGap = 1.1e-1;
        // ReplaceOneSide 允许尝试单边替换的最大候选边间距离。
        double replaceOneSideMaxGap = 1.0e-1;
        // SurfaceIntersection 允许处理的最大边间距离。
        double surfaceIntersectionMaxGap = 1.0e-1;
        // SurfaceIntersection 的求交容差；crossed-face 后置修复不使用此参数。
        double surfaceIntersectionFuzzyTol = 1.0e-1;
        // 中线塌缩允许处理的最大边间距离；不参与 pcurve/face/ShapeFix 容差计算。
        double midlineMaxGap = 1.1e-1;
        // crossed-face 使用原支撑面重建时允许的最大候选边间距离。
        // 实际门槛还会至少覆盖两个候选 face 的已有 tolerance。
        double crossedFaceOriginalSurfaceMaxGap = 1.5e-1;
        // 普通策略的 pcurve 修复和校验容差下限；不降低原候选 edge 的已有 tolerance。
        double pcurveTol = 1.0e-1;
        // face 重建容差，用于新建 face 和局部 wire/face 修复。
        double faceTol = 1.0e-5;
        // face 容差上限。
        double maxFaceTol = 1.0e-3;
        // crossed-face 原支撑面重建允许边界 edge 偏离支撑面的最大实测距离。
        double maxBoundaryToSurfaceDeviation = 2.0e-1;
        // crossed-face 原支撑面重建允许边界 edge 达到的最大拓扑容差。
        double maxRebuiltEdgeTolerance = 2.0e-1;
        // ShapeFix_Shape 的工作精度。
        double shapeFixPrecision = 1.0e-6;
        // ShapeFix_Shape 允许的最大容差。
        double shapeFixMaxTolerance = 1.0e-3;
        // 端点几何距离兜底匹配容差；优先使用 TopoDS_Vertex::IsSame。
        double maxEndpointTolerance = 1.0e-7;
        // 候选边切向允许的最大夹角，方向相反也视为可接受。
        double maxTangentAngle = 20.0 * 3.14159265358979323846 / 180.0;
        // 两条边长度相似度下限：min(lengthA, lengthB) / max(lengthA, lengthB)。
        double minLengthRatio = 0.9;
        // 固定采样点数；关闭自适应采样或自适应参数无效时作为兜底值。
        int sampleCount = 20;
        // 是否启用按边长自适应采样；false 时使用固定 sampleCount。
        bool useAdaptiveSampleCount = true;
        // 自适应采样的最小采样点数，至少包含两个端点。
        int minSampleCount = 8;
        // 自适应采样的最大采样点数，防止长边导致投影计算过慢。
        int maxSampleCount = 80;
        // 目标采样间距，单位与模型单位一致。
        double targetSampleSpacing = 1.0;
        // 为 true 时，在容差不明显更优的情况下优先保留较长边。
        bool preferLongerEdgeToKeep = true;
        // 为 true 时，单次局部替换校验通过后再运行 ShapeFix_Shape。
        bool runShapeFixAfterRepair = true;
        // 预留调试开关；当前 messages_ 始终保留关键过程日志。
        bool verbose = false;
    };

    // 使用指定参数创建修复器。
    explicit HealNarrowGap(const Params& params = Params());

    // 设置输入 shape，并清空上一次修复状态。
    void setInputShape(const TopoDS_Shape& inputShape);

    // 获取最近一次修复后的 shape；没有修改时等于输入 shape。
    TopoDS_Shape getResultShape() const { return resultShape_; }

    // 获取最近一次输入的 shape。
    TopoDS_Shape getInputShape() const { return inputShape_; }

    // 最近一次执行是否发生修改。
    bool isModified() const { return modified_; }

    // 最近一次执行成功修复的候选数量。
    int getRepairedCount() const { return repairedCount_; }

    // 最近一次执行的日志。
    const std::vector<std::string>& getMessages() const { return messages_; }

    // 使用当前 inputShape_ 执行修复；返回成功提交的局部修复次数。
    int perform();

    // 便捷入口：先 setInputShape(inputShape)，再执行 Perform()。
    int perform(const TopoDS_Shape& inputShape);

private:
    // 原支撑面修复的结构化结果；NeedsSurfaceRebuild 当前只记录并拒绝。
    enum class OriginalSurfaceRepairStatus {
        Success,
        NeedsSurfaceRebuild,
        RejectCandidate
    };

    // 单个候选边对及其几何评分数据。
    struct Candidate {
        // 候选边 A 及其所属唯一面。
        TopoDS_Edge edgeA;
        // 候选边 B 及其所属唯一面。
        TopoDS_Edge edgeB;
        // edgeA 所属面。
        TopoDS_Face faceA;
        // edgeB 所属面。
        TopoDS_Face faceB;     
        // edgeA 的 3D 长度。
        double lengthA = 0.0;
        // edgeB 的 3D 长度。
        double lengthB = 0.0;
        // 两条 3D 曲线双向最近距离采样得到的平均间隙。
        double meanGap = 0.0;
        // 两条 3D 曲线双向最近距离采样得到的最大间隙。
        double maxGap = 0.0;
        // 本候选实际使用的采样点数，便于调试自适应采样。
        int usedSampleCount = 0;
        // 最差采样点处的 abs(dot(tangentA, tangentB))。
        double tangentScore = 0.0;
        // 长度相似度，越接近 1 越好。
        double lengthRatio = 0.0;
        // 综合置信度，越大越优先尝试。
        double confidence = 0.0;
    };

    // 普通候选阶段未能处理的边对及其失败原因，供 crossed-face 后置筛选复用。
    struct SkippedCandidate {
        Candidate candidate;
        std::string reason;
    };

    // 单次修复时选择保留哪条边、替换哪条边。
    struct EdgeChoice {
        // 被保留并共享给两个面的边。
        TopoDS_Edge keptEdge;
        // 需要从其 face wire 中替换掉的边。
        TopoDS_Edge edgeToReplace;
        // keptEdge 原本所属的面。
        TopoDS_Face fixedFace;
        // edgeToReplace 所属、需要重建的面。
        TopoDS_Face faceToRebuild;
    };

    // 当前修复器参数。
    Params myParams;

    // 当前输入 shape。
    TopoDS_Shape inputShape_;

    // 最近一次修复后的 shape。
    TopoDS_Shape resultShape_;

    // 最近一次执行是否成功修改。
    bool modified_ = false;

    // 最近一次执行成功修复的候选数量。
    int repairedCount_ = 0;

    // 最近一次执行日志。
    std::vector<std::string> messages_;

    // 从当前 shape 中提取满足拓扑和几何条件的候选边对。
    std::vector<Candidate> extractCandidatePairs(
        const TopoDS_Shape& shape,
        const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
        std::vector<SkippedCandidate>& skippedCandidates,
        std::vector<std::string>& messages) const;

    // 从普通阶段跳过的边对中筛选具有穿插证据的 crossed-face 候选。
    std::vector<Candidate> selectCrossedFaceCandidates(
        const TopoDS_Shape& currentShape,
        const std::vector<SkippedCandidate>& skippedCandidates,
        const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
        std::vector<std::string>& messages) const;

    // 判断两条边是否共享同一对端点，允许方向相反。
    bool haveSameVertexPair(const TopoDS_Edge& edgeA, const TopoDS_Edge& edgeB) const;
    // 按 edge 当前方向获取起点和终点。
    bool getEdgeVerticesOriented(const TopoDS_Edge& edge, TopoDS_Vertex& start, TopoDS_Vertex& end) const;
    // 判断是否为退化边。
    bool isDegeneratedEdge(const TopoDS_Edge& edge) const;
    // 判断是否为 seam edge。
    bool isSeamEdge(
        const TopoDS_Edge& edge,
        const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap) const;
    // 判断 edge 是否存在 3D curve。
    bool has3DCurve(const TopoDS_Edge& edge) const;

    // 根据候选边长度计算实际采样点数；关闭自适应或参数无效时返回固定 sampleCount。
    int computeSampleCountForCandidate(const Candidate& candidate) const;

    // 对候选边对做双向最近距离采样，计算 meanGap/maxGap/tangentScore。
    bool evaluateGeometricSimilarity(
        Candidate& candidate,
        double maximumGap,
        const char* maximumGapName,
        bool requireTangentCompatibility,
        std::string& errorMessage) const;

    // 根据间隙、长度、切向和容差计算综合置信度。
    double computeCandidateConfidence(const Candidate& candidate) const;
    // 根据长度和容差策略决定保留哪条边。
    EdgeChoice chooseEdgeToKeep(const Candidate& candidate) const;
    // 检查 face 的唯一外边界能否按真实连通顺序完整遍历。
    bool hasTraversableCrossedFaceBoundary(
        const TopoDS_Face& face,
        const TopoDS_Edge& candidateEdge,
        std::string& errorMessage,
        bool* hasIncompleteTraversal = nullptr) const;
    // 检查两条候选 edge 是否在公共端点之外发生几何相交。
    bool haveCandidateEdgesIntersectionAwayFromEndpoints(
        const Candidate& candidate,
        std::string& evidenceMessage,
        std::string& errorMessage) const;
    // 只读采样 edge 到 face 未裁剪支撑面的最大投影偏差。
    bool measureEdgeToFaceSupportDeviation(
        const TopoDS_Edge& edge,
        const TopoDS_Face& face,
        int sampleCount,
        double& maximumDeviation,
        std::string& errorMessage) const;
    // 优先选择边数较少的 face 重建；边数相同时复用普通保留边规则。
    bool chooseCrossedFaceBoundaryEdgeChoice(
        const Candidate& candidate,
        EdgeChoice& choice,
        std::string& selectionMessage,
        std::string& errorMessage) const;

    // 尝试修复单个候选；失败时不修改 currentShape，由调用方尝试下一个候选。
    bool tryRepairCandidate(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 原有策略：保留一侧 edge，并用它替换另一侧 edge。
    bool tryRepairByReplacingOneSide(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 修复局部穿插 face：当前只使用 defective face 的原支撑面。
    bool tryRepairByRebuildingCrossedFaceBoundary(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 对无法遍历的 defective wire，删除原 face 后在原支撑面上重建干净边界。
    bool tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage,
        bool allowBoundaryVertexMerge = true);

    // 同一 defective face 的多条 canonical boundary edge 必须在同一 trial 中提交。
    bool tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
        const TopoDS_Shape& currentShape,
        const std::vector<Candidate>& candidates,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage,
        bool allowBoundaryVertexMerge);

    // 在独立 trial 中继续使用 defective face 的原支撑面。
    OriginalSurfaceRepairStatus tryRepairCrossedFaceOnOriginalSurface(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        double& maxBoundaryDeviation,
        std::string& errorMessage);

    // 预留的拟合曲面实现；当前 crossed-face 后置筛选不会调用。
    bool tryRepairCrossedFaceWithRebuiltSurface(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& rebuildMethod,
        double& maxBoundaryDeviation,
        std::string& errorMessage);

    // 复制候选所在的共同 shell，并将 candidate 映射到副本，隔离试算修改。
    bool createIsolatedCandidateTrial(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& trialShape,
        Candidate& trialCandidate,
        std::string& errorMessage) const;

    // 复制共同 owner shell，并在同一副本中映射一组 crossed-face 候选。
    bool createIsolatedCandidateGroupTrial(
        const TopoDS_Shape& currentShape,
        const std::vector<Candidate>& candidates,
        TopoDS_Shape& trialShape,
        std::vector<Candidate>& trialCandidates,
        std::string& errorMessage) const;

    // 用 canonical kept edge 替换 defective edge，构造单一、无分支的闭合边界。
    bool buildCorrectedFaceBoundary(
        const EdgeChoice& choice,
        TopoDS_Wire& correctedWire,
        std::vector<TopoDS_Edge>& boundaryEdges,
        std::string& errorMessage) const;

    // 将 edge 投影到 face 支撑面，并按实测偏差设置最小必要 tolerance。
    bool projectEdgeToFaceSurface(
        TopoDS_Edge& edge,
        const TopoDS_Face& face,
        double maximumDeviation,
        double maximumEdgeTolerance,
        double& measuredDeviation,
        double& appliedTolerance,
        std::string& errorMessage) const;

    // 替换 defective face，并验证共享 edge、旧 edge、方向和整体拓扑。
    bool replaceCrossedFaceAndValidate(
        const TopoDS_Shape& trialShape,
        const EdgeChoice& choice,
        const TopoDS_Face& rebuiltFace,
        double pcurveValidationTol,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage) const;

    // 求交线策略：用两个 face 的真实 surface-section edge 替换两侧 gap edge。
    bool tryRepairBySurfaceIntersection(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 从 faceA/faceB 的 section 结果中选择最贴近当前 gap 的交线 edge。
    bool buildSurfaceIntersectionEdgeForCandidate(
        const Candidate& candidate,
        double endpointDistanceTolerance,
        TopoDS_Edge& intersectionEdge,
        std::string& errorMessage) const;

    // 中线策略：由两侧 edge 采样生成中线 edge，并让两个 face 都引用这条新 edge。
    bool tryRepairByMidlineCollapse(
        const TopoDS_Shape& currentShape,
        const RepairDiagnostics::FailureList& existingFailures,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 根据候选两侧 edge 的采样中点构造一条 3D 中线 edge。
    bool buildMidlineEdgeForCandidate(
        const Candidate& candidate,
        TopoDS_Edge& midlineEdge,
        std::string& errorMessage) const;

    // 在 oldFace 的 wire 中查找 oldEdge 的实际 edge use。
    bool findEdgeUseInFaceWire(
        const TopoDS_Face& oldFace,
        const TopoDS_Edge& oldEdge,
        TopoDS_Edge& oldEdgeInWire,
        std::string& errorMessage) const;

    // 重建 faceToRebuild 的边界，用 replacementEdgeUse 替换 edgeToReplace。
    bool rebuildFaceByReplacingEdge(
        const TopoDS_Face& faceToRebuild,
        const TopoDS_Edge& edgeToReplace,
        const TopoDS_Edge& replacementEdgeUse,
        double faceTol,
        TopoDS_Face& rebuiltFace,
        std::string& errorMessage) const;

    // 确保 edge 在 face 的 surface 上有可用 pcurve。
    bool ensurePCurveOnFace(
        TopoDS_Edge& edge,
        const TopoDS_Face& face,
        double tol,
        std::string& errorMessage) const;

    // 校验替换后的整体 shape 有效，且 free edge 数量没有增加。
    bool validateWholeShape(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& repairedShape,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    // 复用本轮 currentShape 的错误基线，只检查修复结果。
    bool validateWholeShape(
        const RepairDiagnostics::FailureList& existingFailures,
        const TopoDS_Shape& repairedShape,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    // 校验 edge 级替换后 keptEdge 已经被多个 face 共享，并且 pcurve 可用。
    bool validateEdgeReplacement(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& repairedShape,
        const TopoDS_Edge& keptEdge,
        int freeEdgesBefore,
        double pcurveValidationTol,
        std::string& errorMessage) const;

    // 普通候选的最终校验复用本轮基线；中间 trial 仍使用 beforeShape 重载。
    bool validateEdgeReplacement(
        const RepairDiagnostics::FailureList& existingFailures,
        const TopoDS_Shape& repairedShape,
        const TopoDS_Edge& keptEdge,
        int freeEdgesBefore,
        double pcurveValidationTol,
        std::string& errorMessage) const;

    // 两个整体校验入口共用原有的共享 edge 和 pcurve 检查。
    bool validateSharedEdgeAfterReplacement(
        const TopoDS_Shape& repairedShape,
        const TopoDS_Edge& keptEdge,
        double pcurveValidationTol,
        std::string& errorMessage) const;

    // 追加修复日志。
    void appendMessage(std::vector<std::string>& messages, const std::string& message) const;
};
