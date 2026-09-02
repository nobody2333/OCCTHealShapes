#pragma once

#include <string>
#include <vector>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <Geom2d_Curve.hxx>

class BRepFill_Filling;

// 修复同一 shell 上由 free-edge 闭环围成的简单孔洞。
// 补面依次尝试平面、B-spline G1、B-spline C0 和 owner shell 局部 sewing。
class HealSimpleHole
{
public:
    struct Params
    {
        // 闭环至少包含的 edge 数；当前非局部边界距离分类至少需要四条边。
        int minBoundaryEdgeCount = 4;

        // 闭环允许包含的最大 edge 数，避免自动补面处理过于复杂的边界。
        int maxBoundaryEdgeCount = 200;

        // 边界至少应接触多少个不同的已有 face，默认跳过单个孤立 face 的外轮廓。
        int minOwnerFaceCount = 2;

        // 包围盒对角线小于该值时跳过；可用于排除应交给细缝算法的小闭环。
        double minBoundaryDiameter = 1.0e-1;

        // 包围盒对角线大于该值时跳过，单位与模型一致；小于等于 0 表示不设上限。
        double maxBoundaryDiameter = 1.0e1;

        // 边界周长范围；小于等于 0 的上限表示不限制。
        double minBoundaryPerimeter = 0.0;
        double maxBoundaryPerimeter = 0.0;

        // 任意相邻 edge 在公共顶点处允许的最小夹角；更小的尖角按闭合细缝处理。
        double minBoundaryCornerAngle = 45.0 * 3.14159265358979323846 / 180.0;

        // 每条边内部用于边界自距离分析的采样点数，不包含 edge 两端点。
        int boundarySampleCountPerEdge = 5;

        // 投影时排除 source edge 附近的 loop 局部范围，避免把局部相邻段当作对侧边界。
        double localEdgeExclusionRatio = 0.2;

        // 使用所有“采样点到非局部边界”的最近距离分布中的该分位数作为孔洞间距。
        double boundarySeparationPercentile = 0.25;

        // 孔洞间距必须同时满足绝对下限和相对于边界直径的比例下限。
        double minBoundarySeparation = 2.0e-1;
        double minBoundarySeparationRatio = 5.0e-2;

        // 补面面积范围；小于等于 0 的上限表示不限制。
        double minPatchArea = 0.0;
        double maxPatchArea = 0.0;

        // 所有边界 edge 到拟合平面的最大允许距离。
        double planarityTolerance = 1.0e-4;

        // 是否按推荐顺序启用各补面策略。
        bool tryPlanarPatch = true;
        bool tryBSplineG1Patch = true;
        bool tryBSplineC0Patch = true;
        bool useLocalSewingFallback = true;

        // BRepFill_Filling 的边界逼近参数。
        double fillingTolerance = 1.0e-3;
        // Plate 转换为最终 B-spline surface 时的内部近似容差。
        // 它与 Plate 构造容差分开，最终孔边仍由实测偏差上限独立约束。
        double surfaceApproximationTolerance = 5.0e-3;
        // 最终 B-spline 与整个 Plate 曲面的最大允许偏差；近似器即使未达到
        // 目标容差也可能返回 surface，因此必须单独检查实际 ApproxError。
        double maxSurfaceApproximationDeviation = 4.0e-2;
        double maxFillingG0Error = 2.0e-2;
        // 原 boundary edge 到最终重建补面的最大实测距离，与 Plate G0 上限分开控制。
        double maxRebuiltBoundaryDeviation = 4.0e-2;
        double maxFillingG1Error = 15.0 * 3.14159265358979323846 / 180.0;
        double fillingAngularTolerance = 5.0 * 3.14159265358979323846 / 180.0;
        double fillingCurvatureTolerance = 1.0e-1;
        int fillingDegree = 3;
        int fillingPointsOnCurve = 15;
        int fillingIterations = 2;
        // B-spline 近似从最低复杂度开始逐级尝试；优先保持低 degree，
        // 同一 degree 下按固定步长增加 segments。
        int approximationMinDegree = 3;
        int approximationMinSegments = 8;
        int approximationDegreeStep = 1;
        int approximationSegmentStep = 8;
        int fillingMaxDegree = 3;
        int fillingMaxSegments = 32;

        // 普通复杂度范围失败后，是否继续按相同步长扩展 C0 搜索范围。
        bool tryHighComplexityC0Fallback = true;
        int fallbackFillingMaxDegree = 7;
        int fallbackFillingMaxSegments = 100;
        // 只给最低 degree 额外的渐进分段预算；仍受控制点硬上限约束。
        int lowestDegreeMaxSegments = 160;

        // 最终 B-spline surface 的复杂度硬上限，防止分段预算未用满时仍产生异常大曲面。
        int maxFillingPolesPerDirection = 70;
        int maxFillingControlPoints = 4900;

        // 仅用于最终兜底的 owner shell 局部 sewing，不参与 filling 曲面拟合。
        double localSewingTolerance = 1.0e-3;
        double localSewingMinTolerance = 1.0e-7;
        double localSewingMaxTolerance = 1.0e-2;

        // 新 face 上 pcurve 投影精度和算法新增容差上限；已有边容差不会被该上限反向拒绝。
        double pcurveTolerance = 1.0e-5;
        double maxPCurveTolerance = 1.0e-3;
        // 优先把原始 3D 孔边重新投影到最终 B-spline surface；同时测量
        // Plate pcurve，逐边选择实际偏差较小者，避免低 degree 被旧参数化误判。
        bool reprojectBoundaryOnFinalSurface = true;
        int finalBoundaryProjectionSamples = 41;

        // 新 face 加入 shell 后是否用 ShapeFix_Shell 统一 face orientation。
        bool fixShellOrientation = true;

        // ShapeFix_Shell 的工作精度和最大容差。
        double shellFixPrecision = 1.0e-6;
        double shellFixMaxTolerance = 1.0e-3;

        // 为 true 时，修复结果不得增加 free edge。
        bool rejectFreeEdgeIncrease = true;

        // 为 true 时，修复结果不得引入新的 BRepCheck 错误类型/数量。
        bool rejectInvalidResult = true;

        // 单次 perform 最多提交的补面数量，防止异常模型无限重扫。
        int maxRepairCount = 50;

        bool verbose = false;
    };

    explicit HealSimpleHole(const Params& params = Params());

    void setInputShape(const TopoDS_Shape& shape);

    TopoDS_Shape getInputShape() const { return inputShape_; }
    TopoDS_Shape getResultShape() const { return resultShape_; }
    bool isModified() const { return modified_; }
    int getRepairedCount() const { return repairedCount_; }
    const std::vector<std::string>& getMessages() const { return messages_; }

    int perform();
    int perform(const TopoDS_Shape& shape);

private:
    struct Candidate
    {
        TopoDS_Wire boundaryWire;
        std::vector<TopoDS_Edge> boundaryEdges;
        std::vector<TopoDS_Vertex> boundaryVertices;
        std::vector<TopoDS_Face> boundaryOwnerFaces;
        std::vector<TopoDS_Face> ownerFaces;
        TopoDS_Shell ownerShell;
        double perimeter = 0.0;
        double diameter = 0.0;
        double minimumCornerAngle = 0.0;
        double minimumSampledSeparation = 0.0;
        double meanSampledSeparation = 0.0;
        double boundarySeparation = 0.0;
        double boundarySeparationRatio = 0.0;
        double confidence = 0.0;
    };

    Params params_;
    TopoDS_Shape inputShape_;
    TopoDS_Shape resultShape_;
    bool modified_ = false;
    int repairedCount_ = 0;
    std::vector<std::string> messages_;

    // 自行从单面使用 edge 中组装拓扑闭环，不使用 ShapeAnalysis_FreeBounds。
    std::vector<Candidate> extractCandidates(const TopoDS_Shape& shape);

    // 根据公共顶点夹角和非局部边界投影距离排除闭合细缝类型的 loop。
    bool evaluateHoleBoundary(
        Candidate& candidate,
        std::string& errorMessage) const;

    // 用候选闭环拟合平面，并复用原有 boundary edge 构造补面。
    bool buildPlanarPatchFace(
        const Candidate& candidate,
        TopoDS_Face& patchFace,
        double& patchArea,
        std::string& errorMessage) const;

    // 使用 BRepFill_Filling 构造临时 B-spline face，并检查 G0/G1 误差。
    bool buildTemporaryBSplinePatchFace(
        const Candidate& candidate,
        bool useG1Continuity,
        BRepFill_Filling& filling,
        TopoDS_Face& temporaryPatchFace,
        std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
        double& g0Error,
        double& g1Error,
        std::string& errorMessage,
        int maxDegree,
        int maxSegments) const;

    // 使用严格轮廓准则近似 C0 Plate，避免 BRepFill_Filling 的 10 倍 G0 放大。
    bool buildStrictC0BSplinePatchFace(
        const Candidate& candidate,
        TopoDS_Face& temporaryPatchFace,
        std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
        double& g0Error,
        double& g1Error,
        std::string& errorMessage,
        int maxDegree,
        int maxSegments) const;

    // 将原 boundary edge 的 3D curve 投影到 filling surface，再用原 wire 重建 face。
    bool rebuildBSplinePatchWithOriginalBoundary(
        const Candidate& candidate,
        const TopoDS_Face& temporaryPatchFace,
        const std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
        double g0Error,
        TopoDS_Face& patchFace,
        double& patchArea,
        int& reprojectedEdgeCount,
        int& reusedPlatePCurveCount,
        double& maxSelectedPCurveDeviation,
        std::string& errorMessage) const;

    // 修复、校验 patch face 的 pcurve、edge 复用关系和面积。
    bool finalizePatchFace(
        const Candidate& candidate,
        TopoDS_Face& patchFace,
        double baseTolerance,
        double& patchArea,
        std::string& errorMessage) const;

    // 重建候选所属 shell，将补面作为新 face 加入，然后替换原 shell。
    bool addPatchFaceToOwnerShell(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        const TopoDS_Face& patchFace,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage) const;

    // 兜底：只 sewing owner shell 和临时 C0 filling face，并跟踪 sewing 后的共享边。
    bool tryRepairByLocalSewing(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 用已经构造完成的 patch face执行 shell 重建及统一验证。
    bool tryCommitPatchFace(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        const TopoDS_Face& patchFace,
        const std::vector<TopoDS_Edge>& resultBoundaryEdges,
        int freeEdgesBefore,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage) const;

    bool tryRepairCandidate(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& repairedShape,
        std::string& errorMessage);

    // 为单次策略复制 candidate 所属 shell，并把候选映射到复制后的独立拓扑。
    // 策略失败时直接丢弃 trial，避免修改 currentShape 的共享 TShape。
    bool createTrialCandidate(
        const TopoDS_Shape& currentShape,
        const Candidate& candidate,
        TopoDS_Shape& trialShape,
        Candidate& trialCandidate,
        std::string& errorMessage) const;

    bool validateRepairResult(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& afterShape,
        const Candidate& candidate,
        const TopoDS_Face& patchFace,
        const std::vector<TopoDS_Edge>& resultBoundaryEdges,
        int freeEdgesBefore,
        std::string& errorMessage) const;

    void appendMessage(const std::string& message);
};
