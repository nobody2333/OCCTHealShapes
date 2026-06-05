#pragma once

#include <string>
#include <vector>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// 专门处理“同端点双边细缝”的局部修复器。
class HealNarrowGap {
public:
    // 修复参数：只约束本类负责的窄缝场景，不做全局 sewing 或补面。
    struct Params {
        // 两条候选边采样点之间允许的最大 3D 间隙。
        double maxGap = 1.0e-1;
        // 端点几何距离兜底匹配容差；优先使用 TopoDS_Vertex::IsSame。
        double maxEndpointTolerance = 1.0e-7;
        // 候选边切向允许的最大夹角，方向相反也视为可接受。
        double maxTangentAngle = 20.0 * 3.14159265358979323846 / 180.0;
        // 两条边长度相似度下限：min(lengthA, lengthB) / max(lengthA, lengthB)。
        double minLengthRatio = 0.9;
        // 沿 3D 曲线同步采样的点数。
        int sampleCount = 20;
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
        // 两条 3D 曲线同步采样得到的平均间隙。
        double meanGap = 0.0;
        // 两条 3D 曲线同步采样得到的最大间隙。
        double maxGap = 0.0;
        // 最差采样点处的 abs(dot(tangentA, tangentB))。
        double tangentScore = 0.0;
        // 长度相似度，越接近 1 越好。
        double lengthRatio = 0.0;
        // 综合置信度，越大越优先尝试。
        double confidence = 0.0;
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

    // 对候选边对做同步采样，计算 meanGap/maxGap/tangentScore。
    bool evaluateGeometricSimilarity(
        Candidate& candidate,
        std::string& errorMessage) const;

    // 根据间隙、长度、切向和容差计算综合置信度。
    double computeCandidateConfidence(const Candidate& candidate) const;
    // 根据长度和容差策略决定保留哪条边。
    EdgeChoice chooseEdgeToKeep(const Candidate& candidate) const;

    // 尝试修复单个候选；失败时不修改 currentShape，由调用方尝试下一个候选。
    bool tryRepairCandidate(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

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
        double tol,
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

    // 校验 edge 级替换后 keptEdge 已经被多个 face 共享，并且 pcurve 可用。
    bool validateEdgeReplacement(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& repairedShape,
        const TopoDS_Edge& keptEdge,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    // 追加修复日志。
    void appendMessage(std::vector<std::string>& messages, const std::string& message) const;
};
