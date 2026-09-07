#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepLib.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTools_ReShape.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_POnCurv.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Geometry.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <Precision.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TColgp_Array1OfPnt.hxx>

namespace {
constexpr int kMaxAutomaticRepairIterations = 1000;

int computeAutomaticRepairIterationLimit(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) {
        return 1;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);
    if (edgeFaceMap.Extent() <= 0) {
        return 1;
    }

    return std::min(edgeFaceMap.Extent(), kMaxAutomaticRepairIterations);
}

// 将评分项限制在 [0, 1]，避免异常参数影响综合置信度。
double Clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

const char* repairStrategyName(HealNarrowGap::RepairStrategy strategy)
{
    switch (strategy) {
    case HealNarrowGap::RepairStrategy::ReplaceOneSide:
        return "ReplaceOneSide";
    case HealNarrowGap::RepairStrategy::SurfaceIntersection:
        return "SurfaceIntersection";
    case HealNarrowGap::RepairStrategy::MidlineCollapse:
        return "MidlineCollapse";
    case HealNarrowGap::RepairStrategy::Auto:
        return "Auto";
    default:
        return "Unknown";
    }
}

Handle(Geom_Curve) curveForEdgeInWorld(const TopoDS_Edge& edge, double& first, double& last)
{
    first = 0.0;
    last = 0.0;
    if (edge.IsNull()) {
        return Handle(Geom_Curve)();
    }

    TopLoc_Location location;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, location, first, last);
    if (curve.IsNull()) {
        return curve;
    }

    if (!location.IsIdentity()) {
        Handle(Geom_Geometry) transformed = curve->Transformed(location.Transformation());
        curve = Handle(Geom_Curve)::DownCast(transformed);
    }

    return curve;
}

bool projectPointOnCurve(
    const gp_Pnt& point,
    const Handle(Geom_Curve)& curve,
    double first,
    double last,
    double& parameter,
    double& distance)
{
    parameter = first;
    distance = std::numeric_limits<double>::max();
    if (curve.IsNull() || !std::isfinite(first) || !std::isfinite(last) ||
        std::abs(last - first) <= Precision::PConfusion()) {
        return false;
    }

    GeomAPI_ProjectPointOnCurve projector(point, curve, first, last);
    if (projector.NbPoints() <= 0) {
        return false;
    }

    parameter = projector.LowerDistanceParameter();
    distance = projector.LowerDistance();
    return std::isfinite(parameter) && std::isfinite(distance);
}

bool nearestDistanceToEdge(
    const gp_Pnt& point,
    BRepAdaptor_Curve& targetCurve,
    double& nearestDistance)
{
    nearestDistance = std::numeric_limits<double>::max();

    Extrema_ExtPC extrema(
        point,
        targetCurve,
        targetCurve.FirstParameter(),
        targetCurve.LastParameter(),
        Precision::PConfusion());

    if (extrema.IsDone()) {
        for (int extIndex = 1; extIndex <= extrema.NbExt(); ++extIndex) {
            const double squareDistance = extrema.SquareDistance(extIndex);
            if (squareDistance < nearestDistance * nearestDistance) {
                nearestDistance = std::sqrt(std::max(0.0, squareDistance));
            }
        }

        Standard_Real firstSquareDistance = 0.0;
        Standard_Real lastSquareDistance = 0.0;
        gp_Pnt firstPoint;
        gp_Pnt lastPoint;
        extrema.TrimmedSquareDistances(firstSquareDistance, lastSquareDistance, firstPoint, lastPoint);
        nearestDistance = std::min(nearestDistance, std::sqrt(std::max(0.0, firstSquareDistance)));
        nearestDistance = std::min(nearestDistance, std::sqrt(std::max(0.0, lastSquareDistance)));
    }

    if (!std::isfinite(nearestDistance) || nearestDistance == std::numeric_limits<double>::max()) {
        constexpr int fallbackSamples = 32;
        for (int sample = 0; sample < fallbackSamples; ++sample) {
            const double t = static_cast<double>(sample) / static_cast<double>(fallbackSamples - 1);
            const double param = targetCurve.FirstParameter()
                + (targetCurve.LastParameter() - targetCurve.FirstParameter()) * t;
            const double distance = point.Distance(targetCurve.Value(param));
            nearestDistance = std::min(nearestDistance, distance);
        }
    }

    return std::isfinite(nearestDistance) && nearestDistance < std::numeric_limits<double>::max();
}

// 将 OCCT 异常消息转为 std::string，便于写入修复日志。
std::string failureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
}

bool findFacesForSameEdge(const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
                          const TopoDS_Edge& edge,
                          TopTools_ListOfShape& faces)
{
    faces.Clear();
    TopTools_MapOfShape uniqueFaces;
    for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
        const TopoDS_Edge mappedEdge = TopoDS::Edge(edgeFaceMap.FindKey(i));
        if (mappedEdge.IsNull() || !mappedEdge.IsSame(edge)) {
            continue;
        }

        const TopTools_ListOfShape& mappedFaces = edgeFaceMap.FindFromIndex(i);
        for (TopTools_ListIteratorOfListOfShape it(mappedFaces); it.More(); it.Next()) {
            if (uniqueFaces.Add(it.Value())) {
                faces.Append(it.Value());
            }
        }
    }

    return !faces.IsEmpty();
}

bool containsSameEdge(
    const std::vector<TopoDS_Edge>& edges,
    const TopoDS_Edge& target)
{
    for (const TopoDS_Edge& edge : edges) {
        if (!edge.IsNull() && edge.IsSame(target)) {
            return true;
        }
    }
    return false;
}

// 返回 edge 在有向 face 中累计 face orientation 后的唯一 edge-use 方向。
bool findUniqueEdgeUseOrientation(
    const TopoDS_Face& face,
    const TopoDS_Edge& edge,
    TopAbs_Orientation& orientation,
    std::string& errorMessage)
{
    if (face.IsNull() || edge.IsNull()) {
        errorMessage = "face or edge is null while reading edge-use orientation";
        return false;
    }

    int useCount = 0;
    for (TopExp_Explorer edgeExplorer(face, TopAbs_EDGE);
         edgeExplorer.More();
         edgeExplorer.Next()) {
        const TopoDS_Edge edgeUse = TopoDS::Edge(edgeExplorer.Current());
        if (!edgeUse.IsSame(edge)) {
            continue;
        }
        if (edgeUse.Orientation() != TopAbs_FORWARD &&
            edgeUse.Orientation() != TopAbs_REVERSED) {
            errorMessage = "edge use has INTERNAL or EXTERNAL orientation";
            return false;
        }
        orientation = edgeUse.Orientation();
        ++useCount;
    }

    if (useCount != 1) {
        std::ostringstream message;
        message << "edge does not have exactly one use in face: uses=" << useCount;
        errorMessage = message.str();
        return false;
    }
    return true;
}
}

// 构造时规范化参数，避免采样数和迭代数出现无效值。
HealNarrowGap::HealNarrowGap(const Params& params)
    : myParams(params)
{
    myParams.candidateMaxGap = std::max(myParams.candidateMaxGap, Precision::Confusion());
    if (myParams.repairMaxGap <= 0.0) {
        myParams.repairMaxGap = myParams.candidateMaxGap;
    }
    else {
        myParams.repairMaxGap = std::max(myParams.repairMaxGap, Precision::Confusion());
    }
    if (myParams.surfaceIntersectionMaxGap <= 0.0) {
        myParams.surfaceIntersectionMaxGap = myParams.candidateMaxGap;
    }
    else {
        myParams.surfaceIntersectionMaxGap =
            std::max(myParams.surfaceIntersectionMaxGap, Precision::Confusion());
    }
    if (myParams.surfaceIntersectionFuzzyTol < 0.0) {
        myParams.surfaceIntersectionFuzzyTol = 0.0;
    }
    else if (myParams.surfaceIntersectionFuzzyTol > 0.0) {
        myParams.surfaceIntersectionFuzzyTol =
            std::max(myParams.surfaceIntersectionFuzzyTol, Precision::Confusion());
    }
    if (myParams.midlineMaxGap <= 0.0) {
        myParams.midlineMaxGap = myParams.candidateMaxGap;
    }
    else {
        myParams.midlineMaxGap = std::max(myParams.midlineMaxGap, Precision::Confusion());
    }
    if (!std::isfinite(myParams.crossedFaceOriginalSurfaceMaxGap) ||
        myParams.crossedFaceOriginalSurfaceMaxGap <= Precision::Confusion()) {
        myParams.crossedFaceOriginalSurfaceMaxGap = 1.5e-1;
    }

    myParams.pcurveTol = std::max(myParams.pcurveTol, Precision::Confusion());
    myParams.maxPCurveTol = std::max(myParams.maxPCurveTol, myParams.pcurveTol);
    myParams.faceTol = std::max(myParams.faceTol, Precision::Confusion());
    myParams.maxFaceTol = std::max(myParams.maxFaceTol, myParams.faceTol);
    if (!std::isfinite(myParams.maxBoundaryToSurfaceDeviation) ||
        myParams.maxBoundaryToSurfaceDeviation <= Precision::Confusion()) {
        myParams.maxBoundaryToSurfaceDeviation = 1.0e-3;
    }
    if (!std::isfinite(myParams.maxRebuiltEdgeTolerance) ||
        myParams.maxRebuiltEdgeTolerance <= Precision::Confusion()) {
        myParams.maxRebuiltEdgeTolerance = 1.0e-3;
    }
    myParams.maxRebuiltEdgeTolerance = std::max(
        myParams.maxRebuiltEdgeTolerance,
        myParams.pcurveTol);
    myParams.shapeFixPrecision = std::max(myParams.shapeFixPrecision, Precision::Confusion());
    myParams.shapeFixMaxTolerance = std::max(myParams.shapeFixMaxTolerance, myParams.shapeFixPrecision);
    myParams.sampleCount = std::max(2, myParams.sampleCount);
    myParams.minSampleCount = std::max(2, myParams.minSampleCount);
    myParams.maxSampleCount = std::max(myParams.minSampleCount, myParams.maxSampleCount);
    if (myParams.targetSampleSpacing <= Precision::Confusion()) {
        myParams.useAdaptiveSampleCount = false;
    }
    myParams.minLengthRatio = Clamp01(myParams.minLengthRatio);
}

void HealNarrowGap::setInputShape(const TopoDS_Shape& inputShape)
{
    inputShape_ = inputShape;
    resultShape_ = inputShape;
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
}

// 主流程：每轮扫描候选，只提交一个校验通过的最高优先级修复，然后重新扫描。
int HealNarrowGap::perform()
{
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        appendMessage(messages_, "Input shape is null");
        return 0;
    }

    // currentShape 是迭代中的工作副本，不直接修改 inputShape。
    TopoDS_Shape currentShape = inputShape_;
    const int iterationLimit = computeAutomaticRepairIterationLimit(currentShape);
    {
        std::ostringstream os;
        os << "Automatic repair iteration limit = " << iterationLimit;
        appendMessage(messages_, os.str());
    }

    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Original shape", currentShape);

    for (int iter = 0; iter < iterationLimit; ++iter) {

        RepairDiagnostics::AppendIterationHeader(messages_, iter + 1, iter > 0);
        // edgeFaceMap 记录每条 edge 的相邻 face，用于判定单面边和 free edge。
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopologyUtils::BuildEdgeFaceMap(currentShape, edgeFaceMap);
        const bool runCrossedPostPass =
            myParams.repairStrategy == RepairStrategy::Auto &&
            myParams.enableCrossedFaceBoundaryRepair;

        // candidates 是本轮从当前形状重新提取的候选，避免基于过期拓扑批量修改。
        // extractionSkippedCandidates 保存已形成安全边对、但未通过普通几何门槛的候选。
        std::vector<SkippedCandidate> extractionSkippedCandidates;
        std::vector<Candidate> candidates = extractCandidatePairs(
            currentShape,
            edgeFaceMap,
            extractionSkippedCandidates,
            messages_);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "Candidate count = " << candidates.size();
            appendMessage(messages_, os.str());
        }

        const auto appendCandidateFailure = [this](
            const Candidate& candidate,
            const std::string& reason) {
            std::ostringstream failStream;
            failStream << "Candidate failed:"
                       << "\n  maxGap=" << candidate.maxGap
                       << "\n  meanGap=" << candidate.meanGap
                       << "\n  lengthRatio=" << candidate.lengthRatio
                       << "\n  tangentScore=" << candidate.tangentScore
                       << "\n  sampleCount=" << candidate.usedSampleCount
                       << "\n  reason:"
                       << "\n" << RepairDiagnostics::IndentMultiline(
                           RepairDiagnostics::MakeFailureReasonReadable(reason),
                           "    ");
            appendMessage(messages_, failStream.str());
        };

        const auto submitCandidate = [this, &currentShape](
            const Candidate& candidate,
            const TopoDS_Shape& repairedShape,
            const char* submittedStrategy) {
            currentShape = repairedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++repairedCount_;

            std::ostringstream os;
            os << "Repaired candidate: confidence=" << candidate.confidence
               << ", maxGap=" << candidate.maxGap
               << ", meanGap=" << candidate.meanGap
               << ", lengthRatio=" << candidate.lengthRatio
               << ", tangentScore=" << candidate.tangentScore
               << ", sampleCount=" << candidate.usedSampleCount
               << ", repairStrategy=" << submittedStrategy
               << ", candidateMaxGap=" << myParams.candidateMaxGap
               << ", repairMaxGap=" << myParams.repairMaxGap
               << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
               << ", surfaceIntersectionFuzzyTol=" << myParams.surfaceIntersectionFuzzyTol
               << ", midlineMaxGap=" << myParams.midlineMaxGap
               << ", crossedFaceOriginalSurfaceMaxGap="
               << myParams.crossedFaceOriginalSurfaceMaxGap
               << ", shapeFixPrecision=" << myParams.shapeFixPrecision
               << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
            appendMessage(messages_, os.str());
        };

        // 明确存在坏 wire 的候选在常规策略之前保存独立副本。常规策略失败时
        // 可能修改共享 edge 的 pcurve；后置 crossed-face 扫描不能继承这些副作用。
        struct IsolatedCrossedTrial {
            TopoDS_Shape shape;
            Candidate candidate;
            std::string error;
            bool required = false;
            bool available = false;
        };
        std::vector<IsolatedCrossedTrial> crossedTrials(candidates.size());
        if (runCrossedPostPass) {
            for (std::size_t candidateIndex = 0;
                 candidateIndex < candidates.size();
                 ++candidateIndex) {
                const Candidate& candidate = candidates[candidateIndex];
                std::string traversalErrorA;
                std::string traversalErrorB;
                const bool traversableA = hasTraversableCrossedFaceBoundary(
                    candidate.faceA,
                    candidate.edgeA,
                    traversalErrorA);
                const bool traversableB = hasTraversableCrossedFaceBoundary(
                    candidate.faceB,
                    candidate.edgeB,
                    traversalErrorB);
                if (traversableA == traversableB) {
                    continue;
                }

                IsolatedCrossedTrial& trial = crossedTrials[candidateIndex];
                trial.required = true;
                trial.available = createIsolatedCandidateTrial(
                    currentShape,
                    candidate,
                    trial.shape,
                    trial.candidate,
                    trial.error);
            }
        }

        // 第一阶段只运行原有常规策略，并尝试完本轮的全部候选。
        bool repairedThisIteration = false;
        std::vector<std::string> regularFailureReasons(candidates.size());
        for (std::size_t candidateIndex = 0;
             candidateIndex < candidates.size();
             ++candidateIndex) {
            const Candidate& candidate = candidates[candidateIndex];
            TopoDS_Shape repairedShape;
            std::string errorMessage;
            if (!tryRepairCandidate(
                    currentShape,
                    candidate,
                    repairedShape,
                    errorMessage)) {
                regularFailureReasons[candidateIndex] = errorMessage;
                if (!runCrossedPostPass) {
                    appendCandidateFailure(candidate, errorMessage);
                }
                continue;
            }

            submitCandidate(
                candidate,
                repairedShape,
                repairStrategyName(myParams.repairStrategy));
            repairedThisIteration = true;
            break;
        }

        // 只有本轮所有常规候选均未提交时，才从普通阶段跳过的边对中
        // 筛选 crossed-face 候选；这里不再遍历 shape 重新配对。
        if (!repairedThisIteration && runCrossedPostPass) {
            appendMessage(
                messages_,
                "Regular strategies exhausted; selecting CrossedFaceBoundary candidates from regular skips");

            // 同时收集两类普通阶段未处理边对：提取阶段因普通几何门槛跳过的
            // 候选，以及进入三种普通策略但全部失败的候选。
            std::vector<SkippedCandidate> crossedSourceCandidates =
                extractionSkippedCandidates;
            for (std::size_t candidateIndex = 0;
                 candidateIndex < candidates.size();
                 ++candidateIndex) {
                if (regularFailureReasons[candidateIndex].empty()) {
                    continue;
                }
                crossedSourceCandidates.push_back({
                    candidates[candidateIndex],
                    regularFailureReasons[candidateIndex] });
            }

            {
                std::ostringstream message;
                message << "Crossed-face regular skip pool size = "
                        << crossedSourceCandidates.size();
                appendMessage(messages_, message.str());
            }

            std::vector<Candidate> crossedCandidates =
                selectCrossedFaceCandidates(
                    currentShape,
                    crossedSourceCandidates,
                    edgeFaceMap,
                    messages_);
            std::sort(
                crossedCandidates.begin(),
                crossedCandidates.end(),
                [](const Candidate& lhs, const Candidate& rhs) {
                    if (lhs.maxGap != rhs.maxGap) {
                        return lhs.maxGap < rhs.maxGap;
                    }
                    return lhs.confidence > rhs.confidence;
                });

            {
                std::ostringstream message;
                message << "Crossed-face candidate count = "
                        << crossedCandidates.size();
                appendMessage(messages_, message.str());
            }

            const auto sameCandidatePair = [](
                const Candidate& lhs,
                const Candidate& rhs) {
                return (lhs.edgeA.IsSame(rhs.edgeA) &&
                        lhs.edgeB.IsSame(rhs.edgeB)) ||
                       (lhs.edgeA.IsSame(rhs.edgeB) &&
                        lhs.edgeB.IsSame(rhs.edgeA));
            };

            for (const Candidate& sourceCandidate : crossedCandidates) {
                std::size_t regularCandidateIndex = candidates.size();
                for (std::size_t candidateIndex = 0;
                     candidateIndex < candidates.size();
                     ++candidateIndex) {
                    if (sameCandidatePair(
                            sourceCandidate,
                            candidates[candidateIndex])) {
                        regularCandidateIndex = candidateIndex;
                        break;
                    }
                }

                const bool hasRegularCandidate =
                    regularCandidateIndex < candidates.size();
                const IsolatedCrossedTrial* isolatedTrial = hasRegularCandidate
                    ? &crossedTrials[regularCandidateIndex]
                    : nullptr;
                std::string regularFailureReason = hasRegularCandidate
                    ? regularFailureReasons[regularCandidateIndex]
                    : std::string("candidate was skipped by regular candidate filters");
                if (!hasRegularCandidate) {
                    for (const SkippedCandidate& skipped : crossedSourceCandidates) {
                        if (sameCandidatePair(
                                sourceCandidate,
                                skipped.candidate)) {
                            regularFailureReason = skipped.reason;
                            break;
                        }
                    }
                }

                if (isolatedTrial != nullptr &&
                    isolatedTrial->required &&
                    !isolatedTrial->available) {
                    const std::string crossedError =
                        "failed to preserve isolated crossed-face input: "
                        + isolatedTrial->error;
                    appendCandidateFailure(
                        sourceCandidate,
                        regularFailureReason
                        + "\n"
                        + RepairDiagnostics::StrategyFailureBlock(
                            "CrossedFaceBoundary",
                            crossedError));
                    continue;
                }

                const bool useIsolatedTrial =
                    isolatedTrial != nullptr && isolatedTrial->available;
                const TopoDS_Shape& crossedInputShape = useIsolatedTrial
                    ? isolatedTrial->shape
                    : currentShape;
                const Candidate& crossedCandidate = useIsolatedTrial
                    ? isolatedTrial->candidate
                    : sourceCandidate;
                TopoDS_Shape repairedShape;
                std::string crossedError;
                if (!tryRepairByRebuildingCrossedFaceBoundary(
                        crossedInputShape,
                        crossedCandidate,
                        repairedShape,
                        crossedError)) {
                    appendCandidateFailure(
                        sourceCandidate,
                        regularFailureReason
                        + "\n"
                        + RepairDiagnostics::StrategyFailureBlock(
                            "CrossedFaceBoundary",
                            crossedError));
                    continue;
                }

                std::ostringstream warning;
                warning << "Warning: repaired by CrossedFaceBoundary only after all regular "
                        << "candidates failed:"
                        << "\n" << RepairDiagnostics::IndentMultiline(
                            regularFailureReason,
                            "  ");
                appendMessage(messages_, warning.str());
                submitCandidate(
                    sourceCandidate,
                    repairedShape,
                    "CrossedFaceBoundaryPostPass");
                repairedThisIteration = true;
                break;
            }
        }

        if (!repairedThisIteration) {
            appendMessage(messages_, "Repair stopped: no candidate was submitted in this iteration");
            break;
        }

        if (iter + 1 == iterationLimit) {
            appendMessage(messages_, "Repair stopped: automatic iteration limit reached");
        }
    }

    resultShape_ = currentShape;
    appendMessage(messages_, "");
    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Result shape", resultShape_);

    return repairedCount_;
}

int HealNarrowGap::perform(const TopoDS_Shape& inputShape)
{
    setInputShape(inputShape);
    return perform();
}

// 遍历 edgeFaceMap 中的边对，筛出“同端点、各属一个不同面、几何很近”的候选。
std::vector<HealNarrowGap::Candidate> HealNarrowGap::extractCandidatePairs(
    const TopoDS_Shape& shape,
    const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
    std::vector<SkippedCandidate>& skippedCandidates,
    std::vector<std::string>& messages) const
{
    // candidates 保存所有通过拓扑和几何过滤的候选。
    std::vector<Candidate> candidates;
    skippedCandidates.clear();
    if (shape.IsNull()) {
        return candidates;
    }

    for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
        // edgeA 是外层候选边，必须是单面、非 seam、非退化、有 3D curve。
        const TopoDS_Edge edgeA = TopoDS::Edge(edgeFaceMap.FindKey(i));
        if (edgeA.IsNull()) {
            continue;
        }

        // facesA 是 edgeA 的相邻面列表；本修复只处理每条边各属一个面的场景。
        const TopTools_ListOfShape& facesA = edgeFaceMap.FindFromIndex(i);
        if (facesA.Extent() != 1) {
            continue;
        }

        if (isDegeneratedEdge(edgeA) || isSeamEdge(edgeA, edgeFaceMap)) {
            appendMessage(messages, "Candidate skipped: edge is seam or degenerated");
            continue;
        }

        if (!has3DCurve(edgeA)) {
            appendMessage(messages, "Candidate skipped: edge has no 3D curve");
            continue;
        }

        // faceA 是 edgeA 唯一所属面。
        const TopoDS_Face faceA = TopoDS::Face(facesA.First());
        if (faceA.IsNull()) {
            continue;
        }

        for (int j = i + 1; j <= edgeFaceMap.Extent(); ++j) {
            // edgeB 是与 edgeA 配对检查的另一条边。
            const TopoDS_Edge edgeB = TopoDS::Edge(edgeFaceMap.FindKey(j));
            if (edgeB.IsNull() || edgeA.IsSame(edgeB)) {
                continue;
            }

            if (!haveSameVertexPair(edgeA, edgeB)) {
                continue;
            }

            // facesB 是 edgeB 的相邻面列表，同样要求只有一个面。
            const TopTools_ListOfShape& facesB = edgeFaceMap.FindFromIndex(j);
            if (facesB.Extent() != 1) {
                continue;
            }

            // faceB 必须与 faceA 不是同一个拓扑面。
            const TopoDS_Face faceB = TopoDS::Face(facesB.First());
            if (faceB.IsNull() || faceA.IsSame(faceB)) {
                continue;
            }

            if (isDegeneratedEdge(edgeB) || isSeamEdge(edgeB, edgeFaceMap)) {
                appendMessage(messages, "Candidate skipped: edge is seam or degenerated");
                continue;
            }

            if (!has3DCurve(edgeB)) {
                appendMessage(messages, "Candidate skipped: edge has no 3D curve");
                continue;
            }

            // candidate 保存本边对的拓扑引用和后续几何评分。
            Candidate candidate;
            candidate.edgeA = edgeA;
            candidate.edgeB = edgeB;
            candidate.faceA = faceA;
            candidate.faceB = faceB;
            candidate.lengthA = GeometryUtils::ComputeEdgeLength(edgeA);
            candidate.lengthB = GeometryUtils::ComputeEdgeLength(edgeB);

            // maxLength/minLength 用于计算长度相似度。
            const double maxLength = std::max(candidate.lengthA, candidate.lengthB);
            const double minLength = std::min(candidate.lengthA, candidate.lengthB);
            if (maxLength <= Precision::Confusion()) {
                appendMessage(messages, "Candidate skipped: edge length is too small");
                continue;
            }

            candidate.lengthRatio = minLength / maxLength;
            if (candidate.lengthRatio + Precision::Confusion() < myParams.minLengthRatio) {
                appendMessage(messages, "Candidate skipped: length ratio too small");
                continue;
            }

            // errorMessage 记录几何评估失败的具体原因。
            std::string errorMessage;
            if (!evaluateGeometricSimilarity(
                    candidate,
                    myParams.candidateMaxGap,
                    "candidateMaxGap",
                    true,
                    errorMessage)) {
                // 该边对已通过单面边、公共端点、曲线和长度等安全前置检查。
                // 保存后供 crossed-face 使用自己的证据与距离条件再次筛选。
                skippedCandidates.push_back({ candidate, errorMessage });
                appendMessage(messages, "Candidate skipped: " + errorMessage);
                continue;
            }

            candidate.confidence = computeCandidateConfidence(candidate);
            candidates.push_back(candidate);
        }
    }

    return candidates;
}

// 常规阶段耗尽后，只从其跳过列表中筛选 crossed-face 候选。
// 这里不会重新遍历 shape 组合新的 edge 对。
std::vector<HealNarrowGap::Candidate> HealNarrowGap::selectCrossedFaceCandidates(
    const TopoDS_Shape& currentShape,
    const std::vector<SkippedCandidate>& skippedCandidates,
    const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
    std::vector<std::string>& messages) const
{
    std::vector<Candidate> candidates;
    TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
    TopologyUtils::BuildFaceShellMap(currentShape, faceShellMap);

    const auto describeCandidateWire = [](
        const TopoDS_Face& face,
        const TopoDS_Edge& edge) {
        int wireIndex = 0;
        for (TopExp_Explorer wireExplorer(face, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next(), ++wireIndex) {
            int directEdgeUses = 0;
            bool containsCandidate = false;
            for (TopoDS_Iterator edgeIterator(
                     wireExplorer.Current(),
                     Standard_False,
                     Standard_False);
                 edgeIterator.More();
                 edgeIterator.Next()) {
                if (edgeIterator.Value().ShapeType() != TopAbs_EDGE) {
                    continue;
                }
                ++directEdgeUses;
                if (edgeIterator.Value().IsSame(edge)) {
                    containsCandidate = true;
                }
            }
            if (containsCandidate) {
                std::ostringstream description;
                description << wireIndex << ":" << directEdgeUses;
                return description.str();
            }
        }
        return std::string("not-found");
    };

    const auto findEdgeIndex = [&edgeFaceMap](const TopoDS_Edge& edge) {
        for (int edgeIndex = 1;
             edgeIndex <= edgeFaceMap.Extent();
             ++edgeIndex) {
            if (edgeFaceMap.FindKey(edgeIndex).IsSame(edge)) {
                return edgeIndex;
            }
        }
        return 0;
    };

    for (const SkippedCandidate& skipped : skippedCandidates) {
        Candidate candidate = skipped.candidate;
        if (candidate.edgeA.IsNull() || candidate.edgeB.IsNull() ||
            candidate.edgeA.IsSame(candidate.edgeB) ||
            !edgeFaceMap.Contains(candidate.edgeA) ||
            !edgeFaceMap.Contains(candidate.edgeB)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate is stale or no longer belongs to the current shape");
            continue;
        }

        const TopTools_ListOfShape& facesA =
            edgeFaceMap.FindFromKey(candidate.edgeA);
        const TopTools_ListOfShape& facesB =
            edgeFaceMap.FindFromKey(candidate.edgeB);
        if (facesA.Extent() != 1 || facesB.Extent() != 1 ||
            isDegeneratedEdge(candidate.edgeA) ||
            isDegeneratedEdge(candidate.edgeB) ||
            isSeamEdge(candidate.edgeA, edgeFaceMap) ||
            isSeamEdge(candidate.edgeB, edgeFaceMap) ||
            !has3DCurve(candidate.edgeA) ||
            !has3DCurve(candidate.edgeB)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate no longer satisfies safe edge prerequisites");
            continue;
        }

        candidate.faceA = TopoDS::Face(facesA.First());
        candidate.faceB = TopoDS::Face(facesB.First());
        if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
            candidate.faceA.IsSame(candidate.faceB) ||
            !haveSameVertexPair(candidate.edgeA, candidate.edgeB)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate no longer satisfies face or endpoint prerequisites");
            continue;
        }

        // crossed-face 只允许修复同一个唯一 owner shell 内的两个 face。
        // 在 edge 相交和双向投影之前提前排除跨 shell 候选。
        if (!faceShellMap.Contains(candidate.faceA) ||
            !faceShellMap.Contains(candidate.faceB)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: candidate face has no owner shell");
            continue;
        }
        const TopTools_ListOfShape& shellsA =
            faceShellMap.FindFromKey(candidate.faceA);
        const TopTools_ListOfShape& shellsB =
            faceShellMap.FindFromKey(candidate.faceB);
        if (shellsA.Extent() != 1 || shellsB.Extent() != 1) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: candidate face does not have one unique owner shell");
            continue;
        }
        const TopoDS_Shell ownerShellA = TopoDS::Shell(shellsA.First());
        const TopoDS_Shell ownerShellB = TopoDS::Shell(shellsB.First());
        if (ownerShellA.IsNull() || ownerShellB.IsNull() ||
            !ownerShellA.IsSame(ownerShellB)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: candidate faces belong to different owner shells");
            continue;
        }

        // crossed-face 只接受两条候选 edge 在公共端点之外发生相交的边对。
        // face 自身的 wire 分支或遍历异常不再参与重建侧选择。
        std::string intersectionEvidence;
        std::string intersectionError;
        if (!haveCandidateEdgesIntersectionAwayFromEndpoints(
                candidate,
                intersectionEvidence,
                intersectionError)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: "
                    + (intersectionError.empty()
                        ? intersectionEvidence
                        : intersectionError));
            continue;
        }

        EdgeChoice selectionChoice;
        std::string selectionMessage;
        std::string selectionError;
        if (!chooseCrossedFaceBoundaryEdgeChoice(
                candidate,
                selectionChoice,
                selectionMessage,
                selectionError)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: " + selectionError);
            continue;
        }

        candidate.lengthA = GeometryUtils::ComputeEdgeLength(candidate.edgeA);
        candidate.lengthB = GeometryUtils::ComputeEdgeLength(candidate.edgeB);
        const double maxLength = std::max(candidate.lengthA, candidate.lengthB);
        const double minLength = std::min(candidate.lengthA, candidate.lengthB);
        if (maxLength <= Precision::Confusion()) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: edge length is too small");
            continue;
        }

        candidate.lengthRatio = minLength / maxLength;
        if (candidate.lengthRatio + Precision::Confusion() <
            myParams.minLengthRatio) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected: length ratio is too small");
            continue;
        }

        // 该门槛只决定 crossed-face 候选是否进入修复，不参与 pcurve、
        // face 或 ShapeFix 容差计算。
        const double crossedFaceGapLimit = std::max({
            myParams.crossedFaceOriginalSurfaceMaxGap,
            Precision::Confusion(),
            BRep_Tool::Tolerance(candidate.faceA),
            BRep_Tool::Tolerance(candidate.faceB) });
        std::string similarityError;
        if (!evaluateGeometricSimilarity(
                candidate,
                crossedFaceGapLimit,
                "crossedFaceGapLimit",
                false,
                similarityError)) {
            appendMessage(
                messages,
                "Crossed-face skipped candidate rejected after intersection evidence: "
                + similarityError
                + "; intersection="
                + intersectionEvidence
                + "; selection="
                + selectionMessage);
            continue;
        }

        bool isDuplicate = false;
        for (const Candidate& existing : candidates) {
            if ((existing.edgeA.IsSame(candidate.edgeA) &&
                 existing.edgeB.IsSame(candidate.edgeB)) ||
                (existing.edgeA.IsSame(candidate.edgeB) &&
                 existing.edgeB.IsSame(candidate.edgeA))) {
                isDuplicate = true;
                break;
            }
        }
        if (isDuplicate) {
            continue;
        }

        candidate.confidence = computeCandidateConfidence(candidate);
        candidates.push_back(candidate);

        std::ostringstream acceptedMessage;
        acceptedMessage << "Crossed-face candidate selected from regular skips: edgeIndices="
                        << findEdgeIndex(candidate.edgeA) << "/"
                        << findEdgeIndex(candidate.edgeB)
                        << ", sameOwnerShell=true"
                        << ", wireA="
                        << describeCandidateWire(candidate.faceA, candidate.edgeA)
                        << ", wireB="
                        << describeCandidateWire(candidate.faceB, candidate.edgeB)
                        << ", maxGap=" << candidate.maxGap
                        << ", meanGap=" << candidate.meanGap
                        << ", crossedFaceGapLimit=" << crossedFaceGapLimit
                        << ", sampleCount=" << candidate.usedSampleCount
                        << ", intersection=" << intersectionEvidence
                        << ", selection=" << selectionMessage;
        appendMessage(messages, acceptedMessage.str());
    }

    return candidates;
}

// 判断两条边是否引用同一对 TopoDS_Vertex，顺序可相反。
bool HealNarrowGap::haveSameVertexPair(const TopoDS_Edge& edgeA, const TopoDS_Edge& edgeB) const
{
    TopoDS_Vertex aStart;
    TopoDS_Vertex aEnd;
    TopoDS_Vertex bStart;
    TopoDS_Vertex bEnd;
    if (!getEdgeVerticesOriented(edgeA, aStart, aEnd) || !getEdgeVerticesOriented(edgeB, bStart, bEnd)) {
        return false;
    }

    return (aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) ||
           (aStart.IsSame(bEnd) && aEnd.IsSame(bStart));
}

// 获取 edge 在当前 orientation 下的有向起终点。
bool HealNarrowGap::getEdgeVerticesOriented(const TopoDS_Edge& edge, TopoDS_Vertex& start, TopoDS_Vertex& end) const
{
    start.Nullify();
    end.Nullify();
    if (edge.IsNull()) {
        return false;
    }

    TopExp::Vertices(edge, start, end, Standard_True);
    return !start.IsNull() && !end.IsNull();
}

// 退化边不参与窄缝替换，避免破坏特殊面边界。
bool HealNarrowGap::isDegeneratedEdge(const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return true;

    }

    try {
        return BRep_Tool::Degenerated(edge) == Standard_True;
    }
    catch (Standard_Failure&) {
        return true;
    }
}

// seam edge 在周期面上有特殊语义，本类不处理。
bool HealNarrowGap::isSeamEdge(
    const TopoDS_Edge& edge,
    const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap) const
{
    if (edge.IsNull()) {
        return true;
    }

    try {
        if (BRep_Tool::IsClosed(edge)) {
            return true;
        }

        if (!edgeFaceMap.Contains(edge)) {
            return false;
        }

        const TopTools_ListOfShape& faces = edgeFaceMap.FindFromKey(edge);
        for (TopTools_ListIteratorOfListOfShape it(faces); it.More(); it.Next()) {
            const TopoDS_Face face = TopoDS::Face(it.Value());
            if (!face.IsNull() && BRep_Tool::IsClosed(edge, face)) {
                return true;
            }
        }
    }
    catch (Standard_Failure&) {
        return true;
    }

    return false;
}

// 检查 edge 是否已有 3D curve；没有 3D curve 的候选直接跳过。
bool HealNarrowGap::has3DCurve(const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return false;
    }

    try {
        TopLoc_Location location;
        Standard_Real first = 0.0;
        Standard_Real last = 0.0;
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, location, first, last);
        return !curve.IsNull();
    }
    catch (Standard_Failure&) {
        return false;
    }
}

int HealNarrowGap::computeSampleCountForCandidate(const Candidate& candidate) const
{
    const int fixedSampleCount = std::max(2, myParams.sampleCount);
    if (!myParams.useAdaptiveSampleCount) {
        return fixedSampleCount;
    }

    const double length = std::max(candidate.lengthA, candidate.lengthB);
    if (length <= Precision::Confusion()) {
        return fixedSampleCount;
    }

    if (myParams.targetSampleSpacing <= Precision::Confusion()) {
        return fixedSampleCount;
    }

    const double spacing = std::max(myParams.targetSampleSpacing, Precision::Confusion());
    // +1 用于包含首尾端点。
    const int countByLength = static_cast<int>(std::ceil(length / spacing)) + 1;

    return std::max(
        myParams.minSampleCount,
        std::min(myParams.maxSampleCount, countByLength));
}

// 按共享端点在各自 edge 上的曲线参数建立对应采样，评估间隙和切向一致性。
bool HealNarrowGap::evaluateGeometricSimilarity(
    Candidate& candidate,
    double maximumGap,
    const char* maximumGapName,
    bool requireTangentCompatibility,
    std::string& errorMessage) const
{
    try {
        const double safeMaximumGap = std::max(maximumGap, Precision::Confusion());
        const char* safeMaximumGapName = maximumGapName != nullptr
            ? maximumGapName
            : "maximumGap";
        // aStart/aEnd 是 edgeA 在当前拓扑方向下的参考端点。
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        if (!getEdgeVerticesOriented(candidate.edgeA, aStart, aEnd) ||
            !getEdgeVerticesOriented(candidate.edgeB, bStart, bEnd)) {
            errorMessage = "Failed to get oriented edge vertices";
            return false;
        }

        // bParamStartVertex/bParamEndVertex 是 edgeB 上分别对应 aStart/aEnd 的顶点。
        TopoDS_Vertex bParamStartVertex;
        TopoDS_Vertex bParamEndVertex;
        if (aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) {
            bParamStartVertex = bStart;
            bParamEndVertex = bEnd;
        }
        else if (aStart.IsSame(bEnd) && aEnd.IsSame(bStart)) {
            bParamStartVertex = bEnd;
            bParamEndVertex = bStart;
        }
        else {
            errorMessage = "Vertex pair mismatch";
            return false;
        }

        // 参数由对应顶点直接取得，不再用拓扑 orientation 推断几何参数方向。
        const double uA0 = BRep_Tool::Parameter(aStart, candidate.edgeA);
        const double uA1 = BRep_Tool::Parameter(aEnd, candidate.edgeA);
        const double uB0 = BRep_Tool::Parameter(bParamStartVertex, candidate.edgeB);
        const double uB1 = BRep_Tool::Parameter(bParamEndVertex, candidate.edgeB);

        // curveA/curveB 是 OCCT 的边曲线适配器，用于取点和切向。
        BRepAdaptor_Curve curveA(candidate.edgeA);
        BRepAdaptor_Curve curveB(candidate.edgeB);

        // 实际采样点数至少为 2，包含两端点；可按候选边长自适应增加。
        const int sampleCount = computeSampleCountForCandidate(candidate);
        candidate.usedSampleCount = sampleCount;
        // tangentLimit 是允许夹角对应的最小 abs(dot)。
        const double tangentLimit = std::cos(std::max(0.0, myParams.maxTangentAngle));
        // 参数递减采样时，D1 导数方向需要反转成实际采样行进方向。
        const bool reverseTangentA = uA1 < uA0;
        const bool reverseTangentB = uB1 < uB0;

        const auto nearestDistanceToCurve = [](
            const gp_Pnt& point,
            BRepAdaptor_Curve& targetCurve,
            double& nearestParameter,
            double& nearestDistance) {
            nearestParameter = targetCurve.FirstParameter();
            nearestDistance = std::numeric_limits<double>::max();

            Extrema_ExtPC extrema(
                point,
                targetCurve,
                targetCurve.FirstParameter(),
                targetCurve.LastParameter(),
                Precision::PConfusion());

            if (extrema.IsDone()) {
                for (int extIndex = 1; extIndex <= extrema.NbExt(); ++extIndex) {
                    const double squareDistance = extrema.SquareDistance(extIndex);
                    if (squareDistance < nearestDistance * nearestDistance) {
                        nearestDistance = std::sqrt(std::max(0.0, squareDistance));
                        nearestParameter = extrema.Point(extIndex).Parameter();
                    }
                }

                Standard_Real firstSquareDistance = 0.0;
                Standard_Real lastSquareDistance = 0.0;
                gp_Pnt firstPoint;
                gp_Pnt lastPoint;
                extrema.TrimmedSquareDistances(firstSquareDistance, lastSquareDistance, firstPoint, lastPoint);
                if (firstSquareDistance < nearestDistance * nearestDistance) {
                    nearestDistance = std::sqrt(std::max(0.0, firstSquareDistance));
                    nearestParameter = targetCurve.FirstParameter();
                }
                if (lastSquareDistance < nearestDistance * nearestDistance) {
                    nearestDistance = std::sqrt(std::max(0.0, lastSquareDistance));
                    nearestParameter = targetCurve.LastParameter();
                }
            }

            if (!std::isfinite(nearestDistance) || nearestDistance == std::numeric_limits<double>::max()) {
                constexpr int fallbackSamples = 32;
                for (int sample = 0; sample < fallbackSamples; ++sample) {
                    const double t = static_cast<double>(sample) / static_cast<double>(fallbackSamples - 1);
                    const double param = targetCurve.FirstParameter()
                        + (targetCurve.LastParameter() - targetCurve.FirstParameter()) * t;
                    const gp_Pnt candidatePoint = targetCurve.Value(param);
                    const double distance = point.Distance(candidatePoint);
                    if (distance < nearestDistance) {
                        nearestDistance = distance;
                        nearestParameter = param;
                    }
                }
            }

            return std::isfinite(nearestDistance) && nearestDistance < std::numeric_limits<double>::max();
        };

        const auto updateTangentScore = [](
            gp_Vec sourceTangent,
            gp_Vec targetTangent,
            bool reverseSource,
            bool reverseTarget,
            double& minAbsDot) {
            if (reverseSource) {
                sourceTangent.Reverse();
            }
            if (reverseTarget) {
                targetTangent.Reverse();
            }

            if (sourceTangent.SquareMagnitude() <= Precision::SquareConfusion() ||
                targetTangent.SquareMagnitude() <= Precision::SquareConfusion()) {
                return;
            }

            sourceTangent.Normalize();
            targetTangent.Normalize();
            const double absDot = std::abs(sourceTangent.Dot(targetTangent));
            minAbsDot = std::min(minAbsDot, Clamp01(absDot));
        };

        // 双向最近距离可以避免两条 edge 参数分布不一致时的错误对应。
        double bidirectionalGapSum = 0.0;
        double bidirectionalMaxGap = 0.0;
        double maxGapAToB = 0.0;
        double maxGapBToA = 0.0;
        double projectedTangentScore = 1.0;
        int projectedSampleCount = 0;

        for (int i = 0; i < sampleCount; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(sampleCount - 1);
            const double paramA = uA0 + (uA1 - uA0) * t;

            gp_Pnt pointA;
            gp_Vec tangentA;
            curveA.D1(paramA, pointA, tangentA);

            double nearestParamB = 0.0;
            double nearestDistanceB = 0.0;
            if (!nearestDistanceToCurve(pointA, curveB, nearestParamB, nearestDistanceB)) {
                errorMessage = "Failed to project edgeA sample to edgeB";
                return false;
            }

            gp_Pnt projectedPointB;
            gp_Vec projectedTangentB;
            curveB.D1(nearestParamB, projectedPointB, projectedTangentB);
            updateTangentScore(tangentA, projectedTangentB, reverseTangentA, reverseTangentB, projectedTangentScore);

            maxGapAToB = std::max(maxGapAToB, nearestDistanceB);
            bidirectionalMaxGap = std::max(bidirectionalMaxGap, nearestDistanceB);
            bidirectionalGapSum += nearestDistanceB;
            ++projectedSampleCount;

            // 当前方向已有采样超过门槛时，本候选必然失败，无需继续昂贵投影。
            if (nearestDistanceB > safeMaximumGap) {
                std::ostringstream os;
                os << "Candidate gap too large: direction=edgeAtoB"
                   << ", distance=" << nearestDistanceB
                   << ", sample=" << i
                   << ", sampleCount=" << sampleCount
                   << ", " << safeMaximumGapName << "=" << safeMaximumGap;
                errorMessage = os.str();
                return false;
            }
        }

        for (int i = 0; i < sampleCount; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(sampleCount - 1);
            const double paramB = uB0 + (uB1 - uB0) * t;

            gp_Pnt pointB;
            gp_Vec tangentB;
            curveB.D1(paramB, pointB, tangentB);

            double nearestParamA = 0.0;
            double nearestDistanceA = 0.0;
            if (!nearestDistanceToCurve(pointB, curveA, nearestParamA, nearestDistanceA)) {
                errorMessage = "Failed to project edgeB sample to edgeA";
                return false;
            }

            gp_Pnt projectedPointA;
            gp_Vec projectedTangentA;
            curveA.D1(nearestParamA, projectedPointA, projectedTangentA);
            updateTangentScore(tangentB, projectedTangentA, reverseTangentB, reverseTangentA, projectedTangentScore);

            maxGapBToA = std::max(maxGapBToA, nearestDistanceA);
            bidirectionalMaxGap = std::max(bidirectionalMaxGap, nearestDistanceA);
            bidirectionalGapSum += nearestDistanceA;
            ++projectedSampleCount;


            if (nearestDistanceA > safeMaximumGap) {
                std::ostringstream os;
                os << "Candidate gap too large: direction=edgeBtoA"
                   << ", distance=" << nearestDistanceA
                   << ", sample=" << i
                   << ", sampleCount=" << sampleCount
                   << ", " << safeMaximumGapName << "=" << safeMaximumGap;
                errorMessage = os.str();
                return false;
            }
        }

        candidate.meanGap = projectedSampleCount > 0
            ? bidirectionalGapSum / static_cast<double>(projectedSampleCount)
            : 0.0;
        candidate.maxGap = bidirectionalMaxGap;
        candidate.tangentScore = projectedTangentScore;

        if (candidate.maxGap > safeMaximumGap) {
            std::ostringstream os;
            os << "Candidate gap too large: edgeAtoB=" << maxGapAToB
               << ", edgeBtoA=" << maxGapBToA
               << ", sampleCount=" << sampleCount
               << ", " << safeMaximumGapName << "=" << safeMaximumGap;
            errorMessage = os.str();
            return false;
        }

        if (requireTangentCompatibility &&
            candidate.tangentScore + Precision::Confusion() < tangentLimit) {
            errorMessage = "Tangent mismatch";
            return false;
        }

        return true;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Curve sampling failed: " + failureText(failure);
        return false;
    }

    return true;
}

// 将几何指标转换为可排序的综合置信度。
double HealNarrowGap::computeCandidateConfidence(const Candidate& candidate) const
{
    // gapDenom 只使用候选边距离阈值，避免修复容差影响 edge-edge gap 评分。
    const double gapDenom = std::max(myParams.candidateMaxGap, Precision::Confusion());
    // 间隙越小分数越高。
    const double maxGapScore = 1.0 - Clamp01(candidate.maxGap / gapDenom);
    const double meanGapScore = 1.0 - Clamp01(candidate.meanGap / gapDenom);
    // 长度和切向越接近理想值分数越高。
    const double lengthScore = Clamp01(candidate.lengthRatio);
    const double tangentScore = Clamp01(candidate.tangentScore);

    // tolA/tolB 用于给容差更稳定的候选略微加权。
    const double tolA = candidate.edgeA.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeA);
    const double tolB = candidate.edgeB.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeB);
    const double tolSpread = std::abs(tolA - tolB) / std::max(std::max(tolA, tolB), Precision::Confusion());
    const double toleranceScore = 1.0 - Clamp01(tolSpread * 0.25);

    return 0.35 * maxGapScore +
           0.25 * meanGapScore +
           0.20 * lengthScore +
           0.15 * tangentScore +
           0.05 * toleranceScore;
}

// 选择保留边：明显较小容差优先，否则按参数选择较长边或默认 edgeA。
HealNarrowGap::EdgeChoice HealNarrowGap::chooseEdgeToKeep(const Candidate& candidate) const
{
    // choice 输出 keptEdge、edgeToReplace 以及对应两个 face。
    EdgeChoice choice;
    // keepA 表示是否保留 candidate.edgeA。
    bool keepA = true;

    // tolA/tolB 是两条候选边的拓扑容差。
    const double tolA = BRep_Tool::Tolerance(candidate.edgeA);
    const double tolB = BRep_Tool::Tolerance(candidate.edgeB);
    if (tolA * 2.0 < tolB) {
        keepA = true;
    }
    else if (tolB * 2.0 < tolA) {
        keepA = false;
    }
    else if (myParams.preferLongerEdgeToKeep) {
        keepA = candidate.lengthA >= candidate.lengthB;
    }

    if (keepA) {
        choice.keptEdge = candidate.edgeA;
        choice.edgeToReplace = candidate.edgeB;
        choice.fixedFace = candidate.faceA;
        choice.faceToRebuild = candidate.faceB;
    }
    else {
        choice.keptEdge = candidate.edgeB;
        choice.edgeToReplace = candidate.edgeA;
        choice.fixedFace = candidate.faceB;
        choice.faceToRebuild = candidate.faceA;
    }

    return choice;
}

bool HealNarrowGap::hasTraversableCrossedFaceBoundary(
    const TopoDS_Face& face,
    const TopoDS_Edge& candidateEdge,
    std::string& errorMessage,
    bool* hasIncompleteTraversal) const
{
    errorMessage.clear();
    if (hasIncompleteTraversal != nullptr) {
        *hasIncompleteTraversal = false;
    }
    if (face.IsNull() || candidateEdge.IsNull()) {
        errorMessage = "face or candidate edge is null";
        return false;
    }

    try {
        TopoDS_Face forwardFace = face;
        forwardFace.Orientation(TopAbs_FORWARD);

        std::vector<TopoDS_Wire> wires;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            wires.push_back(TopoDS::Wire(wireExplorer.Current()));
        }
        if (wires.size() != 1) {
            std::ostringstream message;
            message << "wire count=" << wires.size();
            errorMessage = message.str();
            return false;
        }

        const TopoDS_Wire& wire = wires.front();
        std::vector<TopoDS_Edge> directEdgeUses;
        for (TopoDS_Iterator edgeIterator(wire, Standard_False, Standard_False);
             edgeIterator.More();
             edgeIterator.Next()) {
            if (edgeIterator.Value().ShapeType() == TopAbs_EDGE) {
                directEdgeUses.push_back(TopoDS::Edge(edgeIterator.Value()));
            }
        }

        std::vector<TopoDS_Edge> orderedEdgeUses;
        for (BRepTools_WireExplorer edgeExplorer(wire, forwardFace);
             edgeExplorer.More();
             edgeExplorer.Next()) {
            orderedEdgeUses.push_back(edgeExplorer.Current());
        }
        if (directEdgeUses.empty() ||
            orderedEdgeUses.size() != directEdgeUses.size()) {
            if (hasIncompleteTraversal != nullptr) {
                *hasIncompleteTraversal = true;
            }
            std::ostringstream message;
            message << "wire traversal is incomplete: directEdgeUses="
                    << directEdgeUses.size()
                    << ", exploredEdgeUses=" << orderedEdgeUses.size();
            errorMessage = message.str();
            return false;
        }

        std::vector<bool> matchedOrderedUses(orderedEdgeUses.size(), false);
        int candidateUseCount = 0;
        for (const TopoDS_Edge& directUse : directEdgeUses) {
            bool matched = false;
            for (std::size_t orderedIndex = 0;
                 orderedIndex < orderedEdgeUses.size();
                 ++orderedIndex) {
                if (!matchedOrderedUses[orderedIndex] &&
                    directUse.IsSame(orderedEdgeUses[orderedIndex])) {
                    matchedOrderedUses[orderedIndex] = true;
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                if (hasIncompleteTraversal != nullptr) {
                    *hasIncompleteTraversal = true;
                }
                errorMessage = "wire traversal omitted a direct edge use";
                return false;
            }
            if (directUse.IsSame(candidateEdge)) {
                ++candidateUseCount;
            }
        }
        if (candidateUseCount != 1) {
            std::ostringstream message;
            message << "candidate edge uses in wire=" << candidateUseCount;
            errorMessage = message.str();
            return false;
        }

    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while checking crossed-face boundary: "
            + failureText(failure);
        return false;
    }
    return true;
}

bool HealNarrowGap::haveCandidateEdgesIntersectionAwayFromEndpoints(
    const Candidate& candidate,
    std::string& evidenceMessage,
    std::string& errorMessage) const
{
    evidenceMessage.clear();
    errorMessage.clear();
    if (candidate.edgeA.IsNull() || candidate.edgeB.IsNull() ||
        candidate.faceA.IsNull() || candidate.faceB.IsNull()) {
        errorMessage = "candidate edge or face is null";
        return false;
    }

    try {
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        TopExp::Vertices(candidate.edgeA, aStart, aEnd, Standard_True);
        TopExp::Vertices(candidate.edgeB, bStart, bEnd, Standard_True);
        if (aStart.IsNull() || aEnd.IsNull() ||
            bStart.IsNull() || bEnd.IsNull()) {
            errorMessage = "candidate edge has a null endpoint";
            return false;
        }

        const double intersectionTolerance = std::max({
            Precision::Confusion(),
            myParams.maxEndpointTolerance,
            BRep_Tool::Tolerance(candidate.faceA),
            BRep_Tool::Tolerance(candidate.faceB),
            BRep_Tool::Tolerance(candidate.edgeA),
            BRep_Tool::Tolerance(candidate.edgeB),
            BRep_Tool::Tolerance(aStart),
            BRep_Tool::Tolerance(aEnd),
            BRep_Tool::Tolerance(bStart),
            BRep_Tool::Tolerance(bEnd) });
        BRepExtrema_DistShapeShape distanceTool(
            candidate.edgeA,
            candidate.edgeB,
            intersectionTolerance);
        if (!distanceTool.IsDone() || !std::isfinite(distanceTool.Value())) {
            errorMessage = "candidate edge intersection check did not finish";
            return false;
        }
        if (distanceTool.Value() > intersectionTolerance) {
            std::ostringstream message;
            message << "candidate edges do not intersect: distance="
                    << distanceTool.Value()
                    << ", tolerance=" << intersectionTolerance;
            evidenceMessage = message.str();
            return false;
        }

        const gp_Pnt sharedEndpointPoints[2] = {
            BRep_Tool::Pnt(aStart),
            BRep_Tool::Pnt(aEnd)
        };
        for (int solutionIndex = 1;
             solutionIndex <= distanceTool.NbSolution();
             ++solutionIndex) {
            const gp_Pnt& pointOnA = distanceTool.PointOnShape1(solutionIndex);
            const gp_Pnt& pointOnB = distanceTool.PointOnShape2(solutionIndex);
            bool isSharedEndpointContact = false;
            for (const gp_Pnt& endpoint : sharedEndpointPoints) {
                if (pointOnA.Distance(endpoint) <= intersectionTolerance &&
                    pointOnB.Distance(endpoint) <= intersectionTolerance) {
                    isSharedEndpointContact = true;
                    break;
                }
            }
            if (!isSharedEndpointContact) {
                std::ostringstream message;
                message << "candidate edges intersect away from their shared endpoints: distance="
                        << distanceTool.Value()
                        << ", tolerance=" << intersectionTolerance
                        << ", solution=" << solutionIndex;
                evidenceMessage = message.str();
                return true;
            }
        }

        std::ostringstream message;
        message << "candidate edges only meet at shared endpoints: solutions="
                << distanceTool.NbSolution()
                << ", tolerance=" << intersectionTolerance;
        evidenceMessage = message.str();
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while checking candidate edge intersection: "
            + failureText(failure);
        return false;
    }
    return false;
}

bool HealNarrowGap::measureEdgeToFaceSupportDeviation(
    const TopoDS_Edge& edge,
    const TopoDS_Face& face,
    int sampleCount,
    double& maximumDeviation,
    std::string& errorMessage) const
{
    maximumDeviation = 0.0;
    errorMessage.clear();
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "edge or face is null while measuring support deviation";
        return false;
    }

    try {
        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) surface =
            BRep_Tool::Surface(face, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "face has no support surface";
            return false;
        }

        BRepAdaptor_Curve curve(edge);
        const double first = curve.FirstParameter();
        const double last = curve.LastParameter();
        if (!std::isfinite(first) || !std::isfinite(last) ||
            std::abs(last - first) <= Precision::PConfusion()) {
            errorMessage = "edge has an invalid curve range";
            return false;
        }

        const int safeSampleCount = std::max(5, sampleCount);
        const gp_Trsf worldToSurface =
            surfaceLocation.Transformation().Inverted();
        const gp_Trsf surfaceToWorld = surfaceLocation.Transformation();
        for (int sampleIndex = 0;
             sampleIndex < safeSampleCount;
             ++sampleIndex) {
            const double ratio = static_cast<double>(sampleIndex) /
                static_cast<double>(safeSampleCount - 1);
            const double parameter = first + (last - first) * ratio;
            const gp_Pnt worldPoint = curve.Value(parameter);
            const gp_Pnt localPoint = worldPoint.Transformed(worldToSurface);
            GeomAPI_ProjectPointOnSurf projector(localPoint, surface);
            if (projector.NbPoints() <= 0) {
                errorMessage = "point projection onto support surface failed";
                return false;
            }

            const gp_Pnt projectedWorld =
                projector.NearestPoint().Transformed(surfaceToWorld);
            maximumDeviation = std::max(
                maximumDeviation,
                worldPoint.Distance(projectedWorld));
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while measuring support deviation: "
            + failureText(failure);
        return false;
    }
    return std::isfinite(maximumDeviation);
}

bool HealNarrowGap::chooseCrossedFaceBoundaryEdgeChoice(
    const Candidate& candidate,
    EdgeChoice& choice,
    std::string& selectionMessage,
    std::string& errorMessage) const
{
    choice = EdgeChoice();
    selectionMessage.clear();
    errorMessage.clear();

    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "candidate face or edge is null while selecting rebuild face";
        return false;
    }

    try {
        TopTools_IndexedMapOfShape faceAEdges;
        TopTools_IndexedMapOfShape faceBEdges;
        TopExp::MapShapes(candidate.faceA, TopAbs_EDGE, faceAEdges);
        TopExp::MapShapes(candidate.faceB, TopAbs_EDGE, faceBEdges);
        if (faceAEdges.IsEmpty() || faceBEdges.IsEmpty()) {
            errorMessage = "candidate face has no boundary edge";
            return false;
        }

        bool rebuildA = false;
        const char* selectionReason = "fewer boundary edges";
        if (faceAEdges.Extent() < faceBEdges.Extent()) {
            rebuildA = true;
        }
        else if (faceBEdges.Extent() < faceAEdges.Extent()) {
            rebuildA = false;
        }
        else {
            // 边数相同时保持现有的保留边规则，使结果不依赖候选枚举顺序。
            const EdgeChoice fallbackChoice = chooseEdgeToKeep(candidate);
            rebuildA = fallbackChoice.faceToRebuild.IsSame(candidate.faceA);
            selectionReason = "equal edge count; kept-edge fallback";
        }

        if (rebuildA) {
            choice.keptEdge = candidate.edgeB;
            choice.edgeToReplace = candidate.edgeA;
            choice.fixedFace = candidate.faceB;
            choice.faceToRebuild = candidate.faceA;
        }
        else {
            choice.keptEdge = candidate.edgeA;
            choice.edgeToReplace = candidate.edgeB;
            choice.fixedFace = candidate.faceA;
            choice.faceToRebuild = candidate.faceB;
        }

        std::ostringstream message;
        message << "rebuild face selected by edge count: faceEdgesA="
                << faceAEdges.Extent()
                << ", faceEdgesB=" << faceBEdges.Extent()
                << ", rebuildFace=" << (rebuildA ? "A" : "B")
                << ", reason=" << selectionReason;
        selectionMessage = message.str();
        return true;
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while selecting crossed-face rebuild face: "
            + failureText(failure);
        return false;
    }
}

// 根据策略尝试提交单个候选修复；不改变主循环“每轮只提交一个候选”的行为。
bool HealNarrowGap::tryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }

    const bool sideGapAllowed = candidate.maxGap <= myParams.repairMaxGap;
    const bool surfaceIntersectionGapAllowed =
        candidate.maxGap <= myParams.surfaceIntersectionMaxGap;
    const bool midlineGapAllowed = candidate.maxGap <= myParams.midlineMaxGap;

    if (myParams.repairStrategy == RepairStrategy::ReplaceOneSide) {
        return tryRepairByReplacingOneSide(currentShape, candidate, repairedShape, errorMessage);
    }

    if (myParams.repairStrategy == RepairStrategy::SurfaceIntersection) {
        if (!surfaceIntersectionGapAllowed) {
            std::ostringstream os;
            os << "surface intersection gap too large: maxGap=" << candidate.maxGap
               << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
               << ", candidateMaxGap=" << myParams.candidateMaxGap;
            errorMessage = os.str();
            return false;
        }

        return tryRepairBySurfaceIntersection(currentShape, candidate, repairedShape, errorMessage);
    }

    if (myParams.repairStrategy == RepairStrategy::MidlineCollapse) {
        if (!midlineGapAllowed) {
            std::ostringstream os;
            os << "midline gap too large: maxGap=" << candidate.maxGap
               << ", midlineMaxGap=" << myParams.midlineMaxGap
               << ", candidateMaxGap=" << myParams.candidateMaxGap;
            errorMessage = os.str();
            return false;
        }

        return tryRepairByMidlineCollapse(currentShape, candidate, repairedShape, errorMessage);
    }

    // Auto 在这里只运行原有常规策略。CrossedFaceBoundary 由 perform() 在
    // 本轮所有常规候选都无法提交后，从普通跳过列表执行后置筛选。
    std::string sideError;
    if (sideGapAllowed) {
        if (tryRepairByReplacingOneSide(
                currentShape,
                candidate,
                repairedShape,
                sideError)) {
            return true;
        }
    }
    else {
        std::ostringstream os;
        os << "repair gap too large for one-side replacement: maxGap=" << candidate.maxGap
           << ", repairMaxGap=" << myParams.repairMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        sideError = os.str();
    }

    std::string surfaceIntersectionError;
    if (surfaceIntersectionGapAllowed) {
        if (tryRepairBySurfaceIntersection(
                currentShape,
                candidate,
                repairedShape,
                surfaceIntersectionError)) {
            if (!sideError.empty()) {
                std::ostringstream warning;
                warning << "Warning: repaired by surface intersection after one-side replacement failed:"
                        << "\n" << RepairDiagnostics::IndentMultiline(
                            RepairDiagnostics::StrategyFailureBlock(
                                "ReplaceOneSide",
                                sideError),
                            "  ");
                appendMessage(messages_, warning.str());
            }
            return true;
        }
    }
    else {
        std::ostringstream message;
        message << "surface intersection gap too large: maxGap=" << candidate.maxGap
                << ", surfaceIntersectionMaxGap="
                << myParams.surfaceIntersectionMaxGap
                << ", candidateMaxGap=" << myParams.candidateMaxGap;
        surfaceIntersectionError = message.str();
    }

    std::string midlineError;
    if (midlineGapAllowed) {
        if (tryRepairByMidlineCollapse(
                currentShape,
                candidate,
                repairedShape,
                midlineError)) {
            std::ostringstream warning;
            warning << "Warning: repaired by midline collapse after earlier strategies failed:"
                    << "\n" << RepairDiagnostics::IndentMultiline(
                        RepairDiagnostics::StrategyFailureBlock(
                            "ReplaceOneSide",
                            sideError),
                        "  ")
                    << "\n" << RepairDiagnostics::IndentMultiline(
                        RepairDiagnostics::StrategyFailureBlock(
                            "SurfaceIntersection",
                            surfaceIntersectionError),
                        "  ");
            appendMessage(messages_, warning.str());
            return true;
        }
    }
    else {
        std::ostringstream message;
        message << "midline gap too large: maxGap=" << candidate.maxGap
                << ", midlineMaxGap=" << myParams.midlineMaxGap;
        midlineError = message.str();
    }

    {
        std::ostringstream message;
        message << "Auto repair failed:"
                << "\n" << RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::StrategyFailureBlock(
                        "ReplaceOneSide",
                        sideError),
                    "  ")
                << "\n" << RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::StrategyFailureBlock(
                        "SurfaceIntersection",
                        surfaceIntersectionError),
                    "  ")
                << "\n" << RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::StrategyFailureBlock(
                        "MidlineCollapse",
                        midlineError),
                    "  ");
        errorMessage = message.str();
    }
    return false;
}

// 尝试提交单个候选修复；任何校验失败都会返回 false，由调用方回滚并尝试下一候选。
bool HealNarrowGap::tryRepairByReplacingOneSide(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    // repairedShape 只有在整个局部替换和校验通过后才赋值。
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }

    if (candidate.maxGap > myParams.repairMaxGap) {
        std::ostringstream os;
        os << "repair gap too large: maxGap=" << candidate.maxGap
           << ", repairMaxGap=" << myParams.repairMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        errorMessage = os.str();
        return false;
    }

    // freeEdgesBefore 用于保证局部修复不会增加自由边数量。
    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    // choice 决定保留哪条边、重建哪个面。
    const EdgeChoice choice = chooseEdgeToKeep(candidate);
    if (choice.keptEdge.IsNull() || choice.edgeToReplace.IsNull() ||
        choice.fixedFace.IsNull() || choice.faceToRebuild.IsNull()) {
        errorMessage = "Invalid edge choice";
        return false;
    }

    // 新共享 edge 的 pcurve 校验继承两条原候选 edge 中较大的 tolerance。
    // 该容差不参与候选 gap 判断、section 求交或 face/ShapeFix 容差计算。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.keptEdge),
        BRep_Tool::Tolerance(choice.edgeToReplace)
        });

    // faceTol 只用于重建 face 和局部 wire/face 修复，不使用 edge-edge gap。
    const double existingFaceTol = std::max(
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.faceToRebuild));
    const double faceTol = std::min(
        std::max(existingFaceTol, myParams.faceTol),
        myParams.maxFaceTol);

    TopoDS_Edge keptEdgeInFixedWire;
    if (!findEdgeUseInFaceWire(choice.fixedFace,
        choice.keptEdge,
        keptEdgeInFixedWire,
        errorMessage)) {
        errorMessage = "Kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeToReplaceInWire;
    if (!findEdgeUseInFaceWire(choice.faceToRebuild,
        choice.edgeToReplace,
        edgeToReplaceInWire,
        errorMessage)) {
        errorMessage = "Replaced edge not found in rebuilt face wire: " + errorMessage;
        return false;
    }

    if (!ensurePCurveOnFace(keptEdgeInFixedWire, choice.fixedFace, pcurveTol, errorMessage)) {
        errorMessage = "Precheck failed: kept edge has invalid pcurve on fixed face: " + errorMessage;
        return false;
    }

    std::string edgeToReplacePCurveError;
    if (!ensurePCurveOnFace(edgeToReplaceInWire, choice.faceToRebuild, pcurveTol, edgeToReplacePCurveError)) {
        // edgeToReplace 后续会从 faceToRebuild 的 wire 中移除。
        // 因此旧 edge use 的 pcurve 失败只记录为 warning，不阻止后续用 keptEdge 重建目标 face。
        std::ostringstream warning;
        warning << "Warning: precheck found invalid pcurve on edge to be replaced; "
                << "continue one-side replacement:"
                << "\n" << RepairDiagnostics::IndentMultiline(
                    RepairDiagnostics::MakeFailureReasonReadable(edgeToReplacePCurveError), "  ");
        appendMessage(messages_, warning.str());
    }

    TopoDS_Face rebuiltFace;
    if (!rebuildFaceByReplacingEdge(
        choice.faceToRebuild,
        choice.edgeToReplace,
        keptEdgeInFixedWire,
        faceTol,
        rebuiltFace,
        errorMessage)) {
        return false;
    }

    // 对于某些候选，rebuiltFace通不过校验，这里先只做记录
	BRepCheck_Analyzer analyzerFace(rebuiltFace);
    if (!analyzerFace.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt face is invalid before shape replacement and final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFace, rebuiltFace));
	}

    // trialShape 是将目标 face 边界重建后的试运行shape。
    TopoDS_Shape trialShape;
    try {
        BRepTools_ReShape reshaper;
        reshaper.Replace(choice.faceToRebuild, rebuiltFace);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "Replaced shape is null";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Shape replacement failed: " + failureText(failure);
        return false;
    }

    if (!findEdgeUseInFaceWire(choice.fixedFace,
        choice.keptEdge,
        keptEdgeInFixedWire,
        errorMessage)) {
        errorMessage = "Kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    if (!findEdgeUseInFaceWire(rebuiltFace,
        choice.keptEdge,
        edgeToReplaceInWire,
        errorMessage)) {
        errorMessage = "Replaced edge not found in rebuilt face wire: " + errorMessage;
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape keptEdgeFaces;
    if (!findFacesForSameEdge(repairedEdgeFaceMap, choice.keptEdge, keptEdgeFaces)) {
        errorMessage = "Kept edge not found after edge replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(keptEdgeFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge keptEdgeInFace;
        if (!findEdgeUseInFaceWire(face, choice.keptEdge, keptEdgeInFace, pcurveRepairError)) {
            errorMessage = "Failed to find kept edge use after replacement: " + pcurveRepairError;
            return false;
        }

        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                keptEdgeInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "Failed to repair kept edge pcurve after replacement: " + pcurveRepairError;
            return false;
        }
    }

    if (myParams.runShapeFixAfterRepair) {
        try {
            const double fixPrecision = std::max(Precision::Confusion(), myParams.shapeFixPrecision);
            const double fixMaxTol = std::max(fixPrecision, myParams.shapeFixMaxTolerance);
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(fixPrecision);
            fixer.SetMaxTolerance(fixMaxTol);
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix failed: " + failureText(failure);
            return false;
        }

    }

    if (!validateEdgeReplacement(
            currentShape,
            trialShape,
            choice.keptEdge,
            freeEdgesBefore,
            pcurveTol,
            errorMessage)) {
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by one-side replacement: repairStrategy="
           << repairStrategyName(myParams.repairStrategy)
           << ", candidateMaxGap=" << myParams.candidateMaxGap
           << ", repairMaxGap=" << myParams.repairMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", surfaceIntersectionFuzzyTol=" << myParams.surfaceIntersectionFuzzyTol
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidate.maxGap=" << candidate.maxGap
           << ", candidate.meanGap=" << candidate.meanGap
           << ", sampleCount=" << candidate.usedSampleCount
           << ", inheritedPCurveTol=" << pcurveTol
           << ", faceTol=" << faceTol
           << ", shapeFixPrecision=" << myParams.shapeFixPrecision
           << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
        appendMessage(messages_, os.str());
    }

    repairedShape = trialShape;
    return true;
}

bool HealNarrowGap::createIsolatedCandidateTrial(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& trialShape,
    Candidate& trialCandidate,
    std::string& errorMessage) const
{
    std::vector<Candidate> trialCandidates;
    if (!createIsolatedCandidateGroupTrial(
            currentShape,
            std::vector<Candidate>{ candidate },
            trialShape,
            trialCandidates,
            errorMessage)) {
        return false;
    }
    trialCandidate = trialCandidates.front();
    return true;
}

bool HealNarrowGap::createIsolatedCandidateGroupTrial(
    const TopoDS_Shape& currentShape,
    const std::vector<Candidate>& candidates,
    TopoDS_Shape& trialShape,
    std::vector<Candidate>& trialCandidates,
    std::string& errorMessage) const
{
    trialShape.Nullify();
    trialCandidates.clear();
    if (currentShape.IsNull() || candidates.empty()) {
        errorMessage = "cannot create an isolated trial from an empty candidate group";
        return false;
    }

    try {
        TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
        TopologyUtils::BuildFaceShellMap(currentShape, faceShellMap);
        TopoDS_Shell ownerShell;
        for (const Candidate& candidate : candidates) {
            if (candidate.edgeA.IsNull() || candidate.edgeB.IsNull() ||
                candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
                !faceShellMap.Contains(candidate.faceA) ||
                !faceShellMap.Contains(candidate.faceB)) {
                errorMessage = "candidate group contains a null or unowned subshape";
                return false;
            }

            const TopTools_ListOfShape& shellsA =
                faceShellMap.FindFromKey(candidate.faceA);
            const TopTools_ListOfShape& shellsB =
                faceShellMap.FindFromKey(candidate.faceB);
            if (shellsA.Extent() != 1 || shellsB.Extent() != 1) {
                errorMessage = "candidate group faces do not have one unique owner shell";
                return false;
            }
            const TopoDS_Shell shellA = TopoDS::Shell(shellsA.First());
            const TopoDS_Shell shellB = TopoDS::Shell(shellsB.First());
            if (shellA.IsNull() || shellB.IsNull() || !shellA.IsSame(shellB)) {
                errorMessage = "candidate group faces do not belong to the same owner shell";
                return false;
            }
            if (ownerShell.IsNull()) {
                ownerShell = shellA;
            }
            else if (!ownerShell.IsSame(shellA)) {
                errorMessage = "candidate group spans multiple owner shells";
                return false;
            }
        }

        // 整组候选只复制一次共同 owner shell，保证所有 edge/face 映射仍然共享。
        BRepBuilderAPI_Copy copier(
            ownerShell,
            Standard_True,
            Standard_False);
        const TopoDS_Shape copiedOwnerShape = copier.Shape();
        if (copiedOwnerShape.IsNull() ||
            copiedOwnerShape.ShapeType() != TopAbs_SHELL) {
            errorMessage = "isolated group trial did not produce an owner shell";
            return false;
        }

        const auto mapSubShape = [&copier](
            const TopoDS_Shape& source,
            TopAbs_ShapeEnum expectedType,
            TopoDS_Shape& mapped) {
            mapped = copier.ModifiedShape(source);
            if (mapped.IsNull() || mapped.ShapeType() != expectedType) {
                return false;
            }
            mapped.Orientation(source.Orientation());
            return true;
        };

        trialCandidates = candidates;
        for (std::size_t candidateIndex = 0;
             candidateIndex < candidates.size();
             ++candidateIndex) {
            const Candidate& source = candidates[candidateIndex];
            Candidate& mapped = trialCandidates[candidateIndex];
            TopoDS_Shape mappedShape;
            if (!mapSubShape(source.edgeA, TopAbs_EDGE, mappedShape)) {
                errorMessage = "cannot map edgeA to the isolated group trial";
                return false;
            }
            mapped.edgeA = TopoDS::Edge(mappedShape);
            if (!mapSubShape(source.edgeB, TopAbs_EDGE, mappedShape)) {
                errorMessage = "cannot map edgeB to the isolated group trial";
                return false;
            }
            mapped.edgeB = TopoDS::Edge(mappedShape);
            if (!mapSubShape(source.faceA, TopAbs_FACE, mappedShape)) {
                errorMessage = "cannot map faceA to the isolated group trial";
                return false;
            }
            mapped.faceA = TopoDS::Face(mappedShape);
            if (!mapSubShape(source.faceB, TopAbs_FACE, mappedShape)) {
                errorMessage = "cannot map faceB to the isolated group trial";
                return false;
            }
            mapped.faceB = TopoDS::Face(mappedShape);
        }

        TopoDS_Shell sourceKey = ownerShell;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Shell replacementValue = TopoDS::Shell(copiedOwnerShape);
        replacementValue.Orientation(TopAbs_FORWARD);
        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull() ||
            trialShape.ShapeType() != currentShape.ShapeType()) {
            errorMessage = "isolated group owner-shell replacement changed the root shape type";
            trialShape.Nullify();
            return false;
        }

        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopologyUtils::BuildEdgeFaceMap(trialShape, edgeFaceMap);
        for (const Candidate& mapped : trialCandidates) {
            TopTools_ListOfShape edgeAFaces;
            TopTools_ListOfShape edgeBFaces;
            if (!findFacesForSameEdge(edgeFaceMap, mapped.edgeA, edgeAFaces) ||
                !findFacesForSameEdge(edgeFaceMap, mapped.edgeB, edgeBFaces) ||
                edgeAFaces.Extent() != 1 || edgeBFaces.Extent() != 1 ||
                !edgeAFaces.First().IsSame(mapped.faceA) ||
                !edgeBFaces.First().IsSame(mapped.faceB)) {
                errorMessage = "candidate group mapping is stale in the isolated trial";
                trialShape.Nullify();
                return false;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while creating isolated candidate group trial: "
            + failureText(failure);
        trialShape.Nullify();
        trialCandidates.clear();
        return false;
    }
    return true;
}

bool HealNarrowGap::buildCorrectedFaceBoundary(
    const EdgeChoice& choice,
    TopoDS_Wire& correctedWire,
    std::vector<TopoDS_Edge>& boundaryEdges,
    std::string& errorMessage) const
{
    correctedWire.Nullify();
    boundaryEdges.clear();
    if (choice.keptEdge.IsNull() || choice.edgeToReplace.IsNull() ||
        choice.fixedFace.IsNull() || choice.faceToRebuild.IsNull()) {
        errorMessage = "crossed-face edge choice is incomplete";
        return false;
    }
    if (!haveSameVertexPair(choice.keptEdge, choice.edgeToReplace)) {
        errorMessage = "crossed-face edges no longer share the same vertex pair";
        return false;
    }

    try {
        TopoDS_Face forwardFace = choice.faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        std::vector<TopoDS_Wire> faceWires;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            faceWires.push_back(TopoDS::Wire(wireExplorer.Current()));
        }
        if (faceWires.size() != 1) {
            std::ostringstream message;
            message << "defective face has unsupported wire count: wires="
                    << faceWires.size();
            errorMessage = message.str();
            return false;
        }

        const TopoDS_Wire sourceWire = faceWires.front();
        std::vector<TopoDS_Edge> directEdgeUses;
        for (TopExp_Explorer edgeExplorer(sourceWire, TopAbs_EDGE);
             edgeExplorer.More();
             edgeExplorer.Next()) {
            directEdgeUses.push_back(TopoDS::Edge(edgeExplorer.Current()));
        }

        std::vector<TopoDS_Edge> orderedEdgeUses;
        for (BRepTools_WireExplorer edgeExplorer(sourceWire, forwardFace);
             edgeExplorer.More();
             edgeExplorer.Next()) {
            orderedEdgeUses.push_back(edgeExplorer.Current());
        }
        if (orderedEdgeUses.size() != directEdgeUses.size()) {
            std::ostringstream message;
            message << "defective wire cannot be traversed completely: directEdgeUses="
                    << directEdgeUses.size()
                    << ", exploredEdgeUses=" << orderedEdgeUses.size();
            errorMessage = message.str();
            return false;
        }

        std::vector<bool> matchedOrderedUses(orderedEdgeUses.size(), false);
        for (const TopoDS_Edge& directUse : directEdgeUses) {
            bool matched = false;
            for (std::size_t orderedIndex = 0;
                 orderedIndex < orderedEdgeUses.size();
                 ++orderedIndex) {
                if (!matchedOrderedUses[orderedIndex] &&
                    directUse.IsSame(orderedEdgeUses[orderedIndex])) {
                    matchedOrderedUses[orderedIndex] = true;
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                errorMessage = "defective wire traversal omitted a direct edge use";
                return false;
            }
        }

        TopoDS_Edge keptEdgeInFixedWire;
        if (!findEdgeUseInFaceWire(
                choice.fixedFace,
                choice.keptEdge,
                keptEdgeInFixedWire,
                errorMessage)) {
            errorMessage = "canonical kept edge has no unique use in fixed face: "
                + errorMessage;
            return false;
        }

        int replacementCount = 0;
        TopoDS_Edge edgeToReplaceInWire;
        for (const TopoDS_Edge& edgeUse : orderedEdgeUses) {
            if (edgeUse.IsSame(choice.edgeToReplace)) {
                edgeToReplaceInWire = edgeUse;
                ++replacementCount;
            }

            if (isDegeneratedEdge(edgeUse) || BRep_Tool::IsClosed(edgeUse) ||
                BRep_Tool::IsClosed(edgeUse, forwardFace)) {
                errorMessage = "defective wire contains a seam, degenerated or self-closed edge";
                return false;
            }
            if (!has3DCurve(edgeUse) ||
                GeometryUtils::ComputeEdgeLength(edgeUse) <= Precision::Confusion()) {
                errorMessage = "defective wire contains an edge without a stable 3d curve";
                return false;
            }
        }
        if (replacementCount != 1 || edgeToReplaceInWire.IsNull()) {
            std::ostringstream message;
            message << "defective edge does not have exactly one wire use: uses="
                    << replacementCount;
            errorMessage = message.str();
            return false;
        }
        if (containsSameEdge(orderedEdgeUses, choice.keptEdge)) {
            errorMessage = "canonical kept edge already exists in the defective wire";
            return false;
        }

        TopoDS_Vertex oldStart;
        TopoDS_Vertex oldEnd;
        TopoDS_Vertex keptStart;
        TopoDS_Vertex keptEnd;
        if (!getEdgeVerticesOriented(edgeToReplaceInWire, oldStart, oldEnd) ||
            !getEdgeVerticesOriented(keptEdgeInFixedWire, keptStart, keptEnd)) {
            errorMessage = "cannot read oriented vertices for crossed-face replacement";
            return false;
        }

        TopoDS_Edge replacementUse = keptEdgeInFixedWire;
        if (keptStart.IsSame(oldStart) && keptEnd.IsSame(oldEnd)) {
            // 当前 kept edge use 已与 defective wire 的旧 edge use 同向。
        }
        else if (keptStart.IsSame(oldEnd) && keptEnd.IsSame(oldStart)) {
            replacementUse.Reverse();
        }
        else {
            errorMessage = "canonical kept edge cannot preserve defective wire connectivity";
            return false;
        }

        BRepBuilderAPI_MakeWire wireMaker;
        for (const TopoDS_Edge& edgeUse : orderedEdgeUses) {
            const TopoDS_Edge useToAdd = edgeUse.IsSame(choice.edgeToReplace)
                ? replacementUse
                : edgeUse;
            if (containsSameEdge(boundaryEdges, useToAdd)) {
                errorMessage = "corrected boundary would contain a repeated or seam edge use";
                return false;
            }
            wireMaker.Add(useToAdd);
            boundaryEdges.push_back(useToAdd);
        }
        if (!wireMaker.IsDone()) {
            errorMessage = "failed to build the corrected crossed-face wire";
            return false;
        }

        correctedWire = wireMaker.Wire();
        correctedWire.Orientation(sourceWire.Orientation());
        if (!BRep_Tool::IsClosed(correctedWire)) {
            errorMessage = "corrected crossed-face wire is not closed";
            correctedWire.Nullify();
            return false;
        }

        BRepCheck_Analyzer wireAnalyzer(correctedWire, Standard_True);
        if (!wireAnalyzer.IsValid()) {
            errorMessage = "corrected crossed-face wire is invalid: "
                + RepairDiagnostics::DescribeBRepCheckFailures(
                    wireAnalyzer,
                    correctedWire);
            correctedWire.Nullify();
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while constructing corrected face boundary: "
            + failureText(failure);
        correctedWire.Nullify();
        boundaryEdges.clear();
        return false;
    }
    return true;
}

bool HealNarrowGap::projectEdgeToFaceSurface(
    TopoDS_Edge& edge,
    const TopoDS_Face& face,
    double maximumDeviation,
    double maximumEdgeTolerance,
    double& measuredDeviation,
    double& appliedTolerance,
    std::string& errorMessage) const
{
    measuredDeviation = 0.0;
    appliedTolerance = 0.0;
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "edge or face is null during boundary projection";
        return false;
    }

    try {
        const double existingTolerance = std::max(
            BRep_Tool::Tolerance(edge),
            Precision::Confusion());
        const double safeMaximumDeviation = std::max(
            maximumDeviation,
            Precision::Confusion());
        const double safeMaximumEdgeTolerance = std::max(
            maximumEdgeTolerance,
            Precision::Confusion());
        // 已有 tolerance 不能被降低，也不属于本次修复新增的误差预算。
        // 当它已高于配置上限时，只允许维持原值，不能继续增大。
        const double allowedTolerance = std::max(
            existingTolerance,
            safeMaximumEdgeTolerance);

        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) surface =
            BRep_Tool::Surface(face, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "defective face has no support surface";
            return false;
        }

        Standard_Real firstParameter = 0.0;
        Standard_Real lastParameter = 0.0;
        Handle(Geom_Curve) projectionCurve =
            BRep_Tool::Curve(edge, firstParameter, lastParameter);
        if (projectionCurve.IsNull() ||
            !std::isfinite(firstParameter) ||
            !std::isfinite(lastParameter) ||
            lastParameter - firstParameter <= Precision::PConfusion()) {
            errorMessage = "boundary edge has no valid 3d curve range";
            return false;
        }

        // BRep_Tool::Curve(edge, first, last) 已位于世界坐标；投影到带 location
        // 的 surface 前，将曲线变换回 surface 的局部坐标系。
        if (!surfaceLocation.IsIdentity()) {
            const Handle(Geom_Geometry) transformed = projectionCurve->Transformed(
                surfaceLocation.Transformation().Inverted());
            projectionCurve = Handle(Geom_Curve)::DownCast(transformed);
            if (projectionCurve.IsNull()) {
                errorMessage = "failed to transform edge curve into surface coordinates";
                return false;
            }
        }

        const double projectionTolerance = std::min(
            allowedTolerance,
            std::max({
                Precision::Confusion(),
                myParams.pcurveTol,
                existingTolerance }));
        ShapeConstruct_ProjectCurveOnSurface projector;
        projector.Init(surface, projectionTolerance);
        Handle(Geom2d_Curve) projectedPCurve;
        if (!projector.Perform(
                projectionCurve,
                firstParameter,
                lastParameter,
                projectedPCurve) &&
            !projector.PerformByProjLib(
                projectionCurve,
                firstParameter,
                lastParameter,
                projectedPCurve)) {
            errorMessage = "failed to project boundary edge onto support surface";
            return false;
        }
        if (projectedPCurve.IsNull()) {
            errorMessage = "boundary projection returned a null pcurve";
            return false;
        }

        BRep_Builder builder;
        builder.UpdateEdge(
            edge,
            projectedPCurve,
            surface,
            surfaceLocation,
            projectionTolerance);
        builder.Range(
            edge,
            surface,
            surfaceLocation,
            firstParameter,
            lastParameter);

        BRepLib_CheckCurveOnSurface checker(edge, face);
        checker.Perform();
        if (!checker.IsDone() || !std::isfinite(checker.MaxDistance())) {
            errorMessage = "failed to measure boundary-to-surface deviation";
            return false;
        }
        measuredDeviation = checker.MaxDistance();
        if (measuredDeviation > safeMaximumDeviation) {
            std::ostringstream message;
            message << "boundary-to-surface deviation exceeds budget: measured="
                    << measuredDeviation
                    << ", maximum=" << safeMaximumDeviation;
            errorMessage = message.str();
            return false;
        }

        const double requiredTolerance = std::max(
            existingTolerance,
            measuredDeviation + Precision::Confusion());
        if (requiredTolerance > allowedTolerance) {
            std::ostringstream message;
            message << "required edge tolerance exceeds configured budget: required="
                    << requiredTolerance
                    << ", maximum=" << allowedTolerance
                    << ", configuredMaximum="
                    << safeMaximumEdgeTolerance
                    << ", existing=" << existingTolerance;
            errorMessage = message.str();
            return false;
        }
        builder.UpdateEdge(edge, requiredTolerance);

        std::string pcurveError;
        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                edge,
                face,
                requiredTolerance,
                &pcurveError)) {
            errorMessage = "projected pcurve failed validation: " + pcurveError;
            return false;
        }

        BRepLib_CheckCurveOnSurface finalChecker(edge, face);
        finalChecker.Perform();
        if (!finalChecker.IsDone() ||
            !std::isfinite(finalChecker.MaxDistance()) ||
            finalChecker.MaxDistance() > safeMaximumDeviation) {
            errorMessage = "repaired pcurve exceeds boundary-to-surface deviation budget";
            return false;
        }
        measuredDeviation = finalChecker.MaxDistance();
        appliedTolerance = std::max(
            BRep_Tool::Tolerance(edge),
            requiredTolerance);
        if (appliedTolerance > allowedTolerance) {
            errorMessage = "pcurve repair increased edge tolerance beyond its budget";
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while projecting boundary edge: "
            + failureText(failure);
        return false;
    }
    return true;
}

bool HealNarrowGap::replaceCrossedFaceAndValidate(
    const TopoDS_Shape& trialShape,
    const EdgeChoice& choice,
    const TopoDS_Face& rebuiltFace,
    double pcurveValidationTol,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (trialShape.IsNull() || rebuiltFace.IsNull()) {
        errorMessage = "crossed-face trial or rebuilt face is null";
        return false;
    }

    BRepCheck_Analyzer faceAnalyzer(rebuiltFace, Standard_True);
    if (!faceAnalyzer.IsValid()) {
        errorMessage = "rebuilt crossed face is invalid: "
            + RepairDiagnostics::DescribeBRepCheckFailures(
                faceAnalyzer,
                rebuiltFace);
        return false;
    }

    try {
        // source key 和 replacement value 都规范为 FORWARD，让 ReShape 在实际
        // shell face use 上保留 defective face 的原 orientation。
        TopoDS_Face sourceKey = choice.faceToRebuild;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Face replacementValue = rebuiltFace;
        replacementValue.Orientation(TopAbs_FORWARD);
        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        repairedShape = reshaper.Apply(trialShape);
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while replacing crossed face: "
            + failureText(failure);
        repairedShape.Nullify();
        return false;
    }

    if (repairedShape.IsNull() ||
        repairedShape.ShapeType() != trialShape.ShapeType()) {
        errorMessage = "crossed-face replacement changed the root shape type";
        repairedShape.Nullify();
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(trialShape);
    if (!validateEdgeReplacement(
            trialShape,
            repairedShape,
            choice.keptEdge,
            freeEdgesBefore,
            pcurveValidationTol,
            errorMessage)) {
        repairedShape.Nullify();
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(repairedShape, edgeFaceMap);
    TopTools_ListOfShape removedEdgeFaces;
    if (findFacesForSameEdge(
            edgeFaceMap,
            choice.edgeToReplace,
            removedEdgeFaces)) {
        errorMessage = "defective edge is still used after crossed-face replacement";
        repairedShape.Nullify();
        return false;
    }

    TopTools_ListOfShape keptEdgeFaces;
    if (!findFacesForSameEdge(edgeFaceMap, choice.keptEdge, keptEdgeFaces) ||
        keptEdgeFaces.Extent() != 2) {
        errorMessage = "canonical kept edge is not shared by exactly two faces";
        repairedShape.Nullify();
        return false;
    }

    TopAbs_Orientation orientations[2] = {
        TopAbs_FORWARD,
        TopAbs_FORWARD
    };
    int faceIndex = 0;
    for (TopTools_ListIteratorOfListOfShape faceIterator(keptEdgeFaces);
         faceIterator.More();
         faceIterator.Next(), ++faceIndex) {
        std::string orientationError;
        if (!findUniqueEdgeUseOrientation(
                TopoDS::Face(faceIterator.Value()),
                choice.keptEdge,
                orientations[faceIndex],
                orientationError)) {
            errorMessage = "cannot validate canonical edge orientation: "
                + orientationError;
            repairedShape.Nullify();
            return false;
        }
    }
    if (orientations[0] == orientations[1]) {
        errorMessage = "canonical kept edge has the same orientation in both faces";
        repairedShape.Nullify();
        return false;
    }
    return true;
}

bool HealNarrowGap::tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage,
    bool allowBoundaryVertexMerge)
{
    return tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
        currentShape,
        std::vector<Candidate>{ candidate },
        repairedShape,
        errorMessage,
        allowBoundaryVertexMerge);
}

bool HealNarrowGap::tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
    const TopoDS_Shape& currentShape,
    const std::vector<Candidate>& candidates,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage,
    bool allowBoundaryVertexMerge)
{
    repairedShape.Nullify();
    if (candidates.empty()) {
        errorMessage = "DeleteAndRebuildOnOriginalSurface received an empty candidate group";
        return false;
    }

    TopoDS_Shape trialShape;
    std::vector<Candidate> trialCandidates;
    if (!createIsolatedCandidateGroupTrial(
            currentShape,
            candidates,
            trialShape,
            trialCandidates,
            errorMessage)) {
        return false;
    }
    const Candidate& trialCandidate = trialCandidates.front();
    const Candidate& candidate = trialCandidate;

    EdgeChoice choice;
    std::string selectionMessage;
    if (!chooseCrossedFaceBoundaryEdgeChoice(
            trialCandidate,
            choice,
            selectionMessage,
            errorMessage)) {
        return false;
    }
    std::string defectiveWireError;
    if (hasTraversableCrossedFaceBoundary(
            choice.faceToRebuild,
            choice.edgeToReplace,
            defectiveWireError)) {
        errorMessage = "DeleteAndRebuildOnOriginalSurface requires an untraversable defective wire";
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, edgeFaceMap);

    struct GroupedChoice {
        Candidate candidate;
        EdgeChoice choice;
    };
    std::vector<GroupedChoice> groupedChoices;
    const auto addGroupedChoice = [&groupedChoices](
        const Candidate& groupedCandidate,
        const EdgeChoice& groupedChoice) {
        for (const GroupedChoice& existing : groupedChoices) {
            if (existing.choice.edgeToReplace.IsSame(groupedChoice.edgeToReplace) ||
                existing.choice.keptEdge.IsSame(groupedChoice.keptEdge)) {
                return;
            }
        }
        groupedChoices.push_back({ groupedCandidate, groupedChoice });
    };
    addGroupedChoice(trialCandidate, choice);
    for (std::size_t candidateIndex = 1;
         candidateIndex < trialCandidates.size();
         ++candidateIndex) {
        const Candidate& groupedCandidate = trialCandidates[candidateIndex];
        EdgeChoice groupedChoice;
        if (groupedCandidate.faceA.IsSame(choice.faceToRebuild) &&
            !groupedCandidate.faceB.IsSame(choice.faceToRebuild)) {
            groupedChoice.keptEdge = groupedCandidate.edgeB;
            groupedChoice.edgeToReplace = groupedCandidate.edgeA;
            groupedChoice.fixedFace = groupedCandidate.faceB;
            groupedChoice.faceToRebuild = groupedCandidate.faceA;
        }
        else if (groupedCandidate.faceB.IsSame(choice.faceToRebuild) &&
                 !groupedCandidate.faceA.IsSame(choice.faceToRebuild)) {
            groupedChoice.keptEdge = groupedCandidate.edgeA;
            groupedChoice.edgeToReplace = groupedCandidate.edgeB;
            groupedChoice.fixedFace = groupedCandidate.faceA;
            groupedChoice.faceToRebuild = groupedCandidate.faceB;
        }
        else {
            errorMessage = "candidate group does not share one defective face";
            return false;
        }
        addGroupedChoice(groupedCandidate, groupedChoice);
    }

    // 顶点合并后的补充边对不一定自身发生穿插，但它与 defective edge 已经
    // 共享端点，且属于已由主候选证明存在 crossed-face 缺陷的同一个 face。
    // 因此这里只补充做长度与双向距离检查，不重复要求相交证据。
    const double groupedGapLimit = std::max({
        myParams.crossedFaceOriginalSurfaceMaxGap,
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.faceToRebuild) });
    for (int defectiveEdgeIndex = 1;
         defectiveEdgeIndex <= edgeFaceMap.Extent();
         ++defectiveEdgeIndex) {
        const TopTools_ListOfShape& defectiveOwners =
            edgeFaceMap.FindFromIndex(defectiveEdgeIndex);
        if (defectiveOwners.Extent() != 1 ||
            !defectiveOwners.First().IsSame(choice.faceToRebuild)) {
            continue;
        }
        const TopoDS_Edge defectiveEdge = TopoDS::Edge(
            edgeFaceMap.FindKey(defectiveEdgeIndex));
        if (defectiveEdge.IsNull() || isDegeneratedEdge(defectiveEdge) ||
            !has3DCurve(defectiveEdge)) {
            continue;
        }

        for (int keptEdgeIndex = 1;
             keptEdgeIndex <= edgeFaceMap.Extent();
             ++keptEdgeIndex) {
            const TopTools_ListOfShape& keptOwners =
                edgeFaceMap.FindFromIndex(keptEdgeIndex);
            if (keptOwners.Extent() != 1 ||
                keptOwners.First().IsSame(choice.faceToRebuild)) {
                continue;
            }
            const TopoDS_Edge keptEdge = TopoDS::Edge(
                edgeFaceMap.FindKey(keptEdgeIndex));
            if (keptEdge.IsNull() || isDegeneratedEdge(keptEdge) ||
                !has3DCurve(keptEdge) ||
                !haveSameVertexPair(defectiveEdge, keptEdge)) {
                continue;
            }

            Candidate supplementalCandidate;
            supplementalCandidate.edgeA = defectiveEdge;
            supplementalCandidate.edgeB = keptEdge;
            supplementalCandidate.faceA = choice.faceToRebuild;
            supplementalCandidate.faceB = TopoDS::Face(keptOwners.First());
            supplementalCandidate.lengthA =
                GeometryUtils::ComputeEdgeLength(defectiveEdge);
            supplementalCandidate.lengthB =
                GeometryUtils::ComputeEdgeLength(keptEdge);
            const double maximumLength = std::max(
                supplementalCandidate.lengthA,
                supplementalCandidate.lengthB);
            if (maximumLength <= Precision::Confusion()) {
                continue;
            }
            supplementalCandidate.lengthRatio = std::min(
                supplementalCandidate.lengthA,
                supplementalCandidate.lengthB) / maximumLength;
            if (supplementalCandidate.lengthRatio + Precision::Confusion() <
                myParams.minLengthRatio) {
                continue;
            }

            std::string supplementalError;
            if (!evaluateGeometricSimilarity(
                    supplementalCandidate,
                    groupedGapLimit,
                    "crossedFaceGroupedGapLimit",
                    false,
                    supplementalError)) {
                continue;
            }
            supplementalCandidate.confidence =
                computeCandidateConfidence(supplementalCandidate);

            EdgeChoice supplementalChoice;
            supplementalChoice.keptEdge = keptEdge;
            supplementalChoice.edgeToReplace = defectiveEdge;
            supplementalChoice.fixedFace = supplementalCandidate.faceB;
            supplementalChoice.faceToRebuild = choice.faceToRebuild;
            addGroupedChoice(supplementalCandidate, supplementalChoice);
        }
    }

    double pcurveTolerance = std::max({
        Precision::Confusion(),
        myParams.pcurveTol,
        BRep_Tool::Tolerance(choice.faceToRebuild) });
    for (const GroupedChoice& grouped : groupedChoices) {
        pcurveTolerance = std::max({
            pcurveTolerance,
            BRep_Tool::Tolerance(grouped.choice.keptEdge),
            BRep_Tool::Tolerance(grouped.choice.edgeToReplace),
            BRep_Tool::Tolerance(grouped.choice.fixedFace) });

        std::string fixedPCurveError;
        if (GeometryUtils::HasValidPCurveOnFace(
                grouped.choice.keptEdge,
                grouped.choice.fixedFace,
                pcurveTolerance,
                &fixedPCurveError)) {
            continue;
        }

        TopoDS_Edge keptEdgeInFixedFace;
        if (!findEdgeUseInFaceWire(
                grouped.choice.fixedFace,
                grouped.choice.keptEdge,
                keptEdgeInFixedFace,
                errorMessage)) {
            errorMessage = "canonical kept edge cannot be located for pcurve repair: "
                + errorMessage;
            return false;
        }

        double fixedBoundaryDeviation = 0.0;
        if (!measureEdgeToFaceSupportDeviation(
                keptEdgeInFixedFace,
                grouped.choice.fixedFace,
                std::max(8, grouped.candidate.usedSampleCount),
                fixedBoundaryDeviation,
                errorMessage)) {
            errorMessage = "canonical kept edge support check failed: " + errorMessage;
            return false;
        }

        const double fixedExistingTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(keptEdgeInFixedFace),
            BRep_Tool::Tolerance(grouped.choice.fixedFace) });
        // 投影检查与输入 tolerance 可能存在极小的数值差异。这里只给已有
        // 精度预算增加 1% 的计算余量，edge tolerance 仍按实测偏差设置。
        const double fixedValidationTolerance = fixedExistingTolerance
            + std::max(Precision::Confusion(), fixedExistingTolerance * 1.0e-2);
        if (fixedBoundaryDeviation > fixedValidationTolerance) {
            std::ostringstream message;
            message << "canonical kept edge exceeds its existing fixed-face tolerance: deviation="
                    << fixedBoundaryDeviation
                    << ", existingTolerance=" << fixedExistingTolerance
                    << ", validationTolerance=" << fixedValidationTolerance;
            errorMessage = message.str();
            return false;
        }

        double fixedAppliedTolerance = 0.0;
        if (!projectEdgeToFaceSurface(
                keptEdgeInFixedFace,
                grouped.choice.fixedFace,
                fixedValidationTolerance,
                fixedValidationTolerance,
                fixedBoundaryDeviation,
                fixedAppliedTolerance,
                errorMessage)) {
            errorMessage = "canonical kept edge pcurve repair failed: " + errorMessage;
            return false;
        }
        pcurveTolerance = std::max(pcurveTolerance, fixedAppliedTolerance);

        std::ostringstream repairedPCurveMessage;
        repairedPCurveMessage
            << "Crossed-face canonical kept-edge pcurve repaired on fixed face: deviation="
            << fixedBoundaryDeviation
            << ", tolerance=" << fixedAppliedTolerance;
        appendMessage(messages_, repairedPCurveMessage.str());
    }

    {
        std::ostringstream groupedMessage;
        groupedMessage << "Crossed-face defective-face group: candidatePairs="
                       << groupedChoices.size();
        appendMessage(messages_, groupedMessage.str());
    }

    try {
        const auto findSingleSurvivingOwner = [&edgeFaceMap, &choice](
            const TopoDS_Edge& edge,
            TopoDS_Face& owner,
            int& ownerCount,
            std::string& ownerError) {
            owner.Nullify();
            ownerCount = 0;
            if (edge.IsNull() || !edgeFaceMap.Contains(edge)) {
                ownerError = "edge is missing from the edge-face map";
                return false;
            }

            const TopTools_ListOfShape& adjacentFaces = edgeFaceMap.FindFromKey(edge);
            for (TopTools_ListIteratorOfListOfShape faceIterator(adjacentFaces);
                 faceIterator.More();
                 faceIterator.Next()) {
                const TopoDS_Face adjacentFace = TopoDS::Face(faceIterator.Value());
                if (adjacentFace.IsNull() ||
                    adjacentFace.IsSame(choice.faceToRebuild)) {
                    continue;
                }
                if (owner.IsNull() || !owner.IsSame(adjacentFace)) {
                    owner = adjacentFace;
                    ++ownerCount;
                }
            }
            if (ownerCount != 1) {
                std::ostringstream message;
                message << "surviving owner faces=" << ownerCount;
                ownerError = message.str();
                owner.Nullify();
                return false;
            }
            return true;
        };

        // 坏 wire 不能依赖 WireExplorer 排序，只读取其直接保存的 edge use。
        TopoDS_Face forwardDefectiveFace = choice.faceToRebuild;
        forwardDefectiveFace.Orientation(TopAbs_FORWARD);
        std::vector<TopoDS_Wire> defectiveWires;
        for (TopExp_Explorer wireExplorer(forwardDefectiveFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            defectiveWires.push_back(TopoDS::Wire(wireExplorer.Current()));
        }
        if (defectiveWires.empty()) {
            errorMessage = "DeleteAndRebuildOnOriginalSurface found no defective wire";
            return false;
        }

        struct BoundaryUse {
            TopoDS_Edge edge;
            TopoDS_Face ownerFace;
        };
        std::vector<BoundaryUse> unorderedBoundary;
        std::vector<TopoDS_Edge> removedEdges;
        const auto addBoundary = [&unorderedBoundary](
            const TopoDS_Edge& edge,
            const TopoDS_Face& ownerFace) {
            for (const BoundaryUse& existing : unorderedBoundary) {
                if (existing.edge.IsSame(edge)) {
                    return false;
                }
            }
            unorderedBoundary.push_back({ edge, ownerFace });
            return true;
        };
        const auto addRemovedEdge = [&removedEdges](const TopoDS_Edge& edge) {
            for (const TopoDS_Edge& existing : removedEdges) {
                if (existing.IsSame(edge)) {
                    return;
                }
            }
            removedEdges.push_back(edge);
        };

        std::vector<int> replacementUses(groupedChoices.size(), 0);
        for (const TopoDS_Wire& defectiveWire : defectiveWires) {
            for (TopoDS_Iterator edgeIterator(
                     defectiveWire,
                     Standard_False,
                     Standard_False);
                 edgeIterator.More();
                 edgeIterator.Next()) {
                if (edgeIterator.Value().ShapeType() != TopAbs_EDGE) {
                    continue;
                }
                const TopoDS_Edge directEdgeUse = TopoDS::Edge(edgeIterator.Value());
                if (isDegeneratedEdge(directEdgeUse) ||
                    BRep_Tool::IsClosed(directEdgeUse) ||
                    BRep_Tool::IsClosed(directEdgeUse, forwardDefectiveFace)) {
                    errorMessage = "defective face contains a seam, degenerated or self-closed edge";
                    return false;
                }
                if (!has3DCurve(directEdgeUse) ||
                    GeometryUtils::ComputeEdgeLength(directEdgeUse) <= Precision::Confusion()) {
                    errorMessage = "defective face contains an edge without a stable 3d curve";
                    return false;
                }

                bool isReplacement = false;
                for (std::size_t choiceIndex = 0;
                     choiceIndex < groupedChoices.size();
                     ++choiceIndex) {
                    if (directEdgeUse.IsSame(
                            groupedChoices[choiceIndex].choice.edgeToReplace)) {
                        ++replacementUses[choiceIndex];
                        isReplacement = true;
                        break;
                    }
                }
                if (isReplacement) {
                    addRemovedEdge(directEdgeUse);
                    continue;
                }

                // 只保留删除 defective face 后仍由另一 face 使用的边。
                TopoDS_Face survivingOwner;
                int ownerCount = 0;
                std::string ownerError;
                if (!findSingleSurvivingOwner(
                        directEdgeUse,
                        survivingOwner,
                        ownerCount,
                        ownerError)) {
                    if (ownerCount == 0) {
                        addRemovedEdge(directEdgeUse);
                        continue;
                    }
                    errorMessage = "cannot preserve defective-face boundary edge: " + ownerError;
                    return false;
                }
                if (!addBoundary(directEdgeUse, survivingOwner)) {
                    errorMessage = "surviving boundary contains a repeated edge use";
                    return false;
                }
            }
        }
        for (std::size_t choiceIndex = 0;
             choiceIndex < groupedChoices.size();
             ++choiceIndex) {
            if (replacementUses[choiceIndex] != 1) {
                std::ostringstream message;
                message << "defective candidate edge uses="
                        << replacementUses[choiceIndex]
                        << ", candidateIndex=" << choiceIndex;
                errorMessage = message.str();
                return false;
            }

            const EdgeChoice& groupedChoice = groupedChoices[choiceIndex].choice;
            TopoDS_Face keptEdgeOwner;
            int keptEdgeOwnerCount = 0;
            std::string keptOwnerError;
            if (!findSingleSurvivingOwner(
                    groupedChoice.keptEdge,
                    keptEdgeOwner,
                    keptEdgeOwnerCount,
                    keptOwnerError) ||
                !keptEdgeOwner.IsSame(groupedChoice.fixedFace)) {
                errorMessage = "canonical kept edge does not have the expected surviving face";
                if (!keptOwnerError.empty()) {
                    errorMessage += ": " + keptOwnerError;
                }
                return false;
            }
            if (!addBoundary(groupedChoice.keptEdge, keptEdgeOwner)) {
                errorMessage = "canonical kept edge is already present in defective face boundary";
                return false;
            }
        }
        if (unorderedBoundary.size() < 3) {
            errorMessage = "DeleteAndRebuildOnOriginalSurface boundary has fewer than three edges";
            return false;
        }

        struct BoundaryLink {
            TopoDS_Edge edge;
            TopoDS_Face ownerFace;
            TopoDS_Vertex firstVertex;
            TopoDS_Vertex lastVertex;
            int firstVertexIndex = 0;
            int lastVertexIndex = 0;
        };
        TopTools_IndexedMapOfShape boundaryVertices;
        std::vector<BoundaryLink> links;
        links.reserve(unorderedBoundary.size());
        for (const BoundaryUse& boundaryUse : unorderedBoundary) {
            BoundaryLink link;
            link.edge = boundaryUse.edge;
            link.ownerFace = boundaryUse.ownerFace;
            if (!getEdgeVerticesOriented(
                    link.edge,
                    link.firstVertex,
                    link.lastVertex) ||
                link.firstVertex.IsSame(link.lastVertex)) {
                errorMessage = "surviving boundary contains an invalid closed edge";
                return false;
            }
            link.firstVertexIndex = boundaryVertices.Add(link.firstVertex);
            link.lastVertexIndex = boundaryVertices.Add(link.lastVertex);
            links.push_back(link);
        }

        std::vector<std::vector<std::size_t>> incidentLinks(
            static_cast<std::size_t>(boundaryVertices.Extent() + 1));
        for (std::size_t linkIndex = 0; linkIndex < links.size(); ++linkIndex) {
            incidentLinks[static_cast<std::size_t>(links[linkIndex].firstVertexIndex)]
                .push_back(linkIndex);
            incidentLinks[static_cast<std::size_t>(links[linkIndex].lastVertexIndex)]
                .push_back(linkIndex);
        }
        std::vector<int> openVertexIndices;
        int invalidVertexCount = 0;
        std::size_t maximumVertexDegree = 0;
        for (int vertexIndex = 1;
             vertexIndex <= boundaryVertices.Extent();
             ++vertexIndex) {
            const std::size_t degree =
                incidentLinks[static_cast<std::size_t>(vertexIndex)].size();
            maximumVertexDegree = std::max(maximumVertexDegree, degree);
            if (degree != 2) {
                ++invalidVertexCount;
            }
            if (degree == 1) {
                openVertexIndices.push_back(vertexIndex);
            }
        }
        if (invalidVertexCount > 0) {
            double nearestOpenVertexDistance = std::numeric_limits<double>::max();
            for (std::size_t firstIndex = 0;
                 firstIndex < openVertexIndices.size();
                 ++firstIndex) {
                const TopoDS_Vertex firstVertex = TopoDS::Vertex(
                    boundaryVertices.FindKey(openVertexIndices[firstIndex]));
                for (std::size_t secondIndex = firstIndex + 1;
                     secondIndex < openVertexIndices.size();
                     ++secondIndex) {
                    const TopoDS_Vertex secondVertex = TopoDS::Vertex(
                        boundaryVertices.FindKey(openVertexIndices[secondIndex]));
                    nearestOpenVertexDistance = std::min(
                        nearestOpenVertexDistance,
                        BRep_Tool::Pnt(firstVertex).Distance(
                            BRep_Tool::Pnt(secondVertex)));
                }
            }

            double bestMissingEdgeEndpointDistance =
                std::numeric_limits<double>::max();
            int bestMissingEdgeIndex = 0;
            int bestMissingEdgeOwnerCount = 0;
            TopoDS_Edge bestMissingEdge;
            TopoDS_Face bestMissingEdgeOwner;
            double supplementalSupportDeviation = 0.0;
            std::string supplementalFailureReason;
            if (openVertexIndices.size() == 2) {
                const TopoDS_Vertex openFirst = TopoDS::Vertex(
                    boundaryVertices.FindKey(openVertexIndices[0]));
                const TopoDS_Vertex openSecond = TopoDS::Vertex(
                    boundaryVertices.FindKey(openVertexIndices[1]));
                const gp_Pnt openFirstPoint = BRep_Tool::Pnt(openFirst);
                const gp_Pnt openSecondPoint = BRep_Tool::Pnt(openSecond);
                for (int edgeIndex = 1;
                     edgeIndex <= edgeFaceMap.Extent();
                     ++edgeIndex) {
                    const TopoDS_Edge candidateBoundaryEdge = TopoDS::Edge(
                        edgeFaceMap.FindKey(edgeIndex));
                    bool alreadyUsed = false;
                    for (const BoundaryUse& boundaryUse : unorderedBoundary) {
                        if (boundaryUse.edge.IsSame(candidateBoundaryEdge)) {
                            alreadyUsed = true;
                            break;
                        }
                    }
                    if (alreadyUsed ||
                        containsSameEdge(removedEdges, candidateBoundaryEdge)) {
                        continue;
                    }

                    TopoDS_Vertex candidateFirst;
                    TopoDS_Vertex candidateSecond;
                    if (!getEdgeVerticesOriented(
                            candidateBoundaryEdge,
                            candidateFirst,
                            candidateSecond)) {
                        continue;
                    }
                    const gp_Pnt firstPoint = BRep_Tool::Pnt(candidateFirst);
                    const gp_Pnt secondPoint = BRep_Tool::Pnt(candidateSecond);
                    const double directDistance = std::max(
                        openFirstPoint.Distance(firstPoint),
                        openSecondPoint.Distance(secondPoint));
                    const double reversedDistance = std::max(
                        openFirstPoint.Distance(secondPoint),
                        openSecondPoint.Distance(firstPoint));
                    const double endpointDistance = std::min(
                        directDistance,
                        reversedDistance);
                    if (endpointDistance < bestMissingEdgeEndpointDistance) {
                        bestMissingEdgeEndpointDistance = endpointDistance;
                        bestMissingEdgeIndex = edgeIndex;
                        bestMissingEdgeOwnerCount =
                            edgeFaceMap.FindFromIndex(edgeIndex).Extent();
                        bestMissingEdge = candidateBoundaryEdge;
                        bestMissingEdgeOwner.Nullify();
                        if (bestMissingEdgeOwnerCount == 1) {
                            bestMissingEdgeOwner = TopoDS::Face(
                                edgeFaceMap.FindFromIndex(edgeIndex).First());
                        }
                    }
                }
            }

            if (allowBoundaryVertexMerge &&
                openVertexIndices.size() == 2 &&
                bestMissingEdgeOwnerCount == 1 &&
                !bestMissingEdge.IsNull() &&
                !bestMissingEdgeOwner.IsNull()) {
                const TopoDS_Vertex openFirst = TopoDS::Vertex(
                    boundaryVertices.FindKey(openVertexIndices[0]));
                const TopoDS_Vertex openSecond = TopoDS::Vertex(
                    boundaryVertices.FindKey(openVertexIndices[1]));
                TopoDS_Vertex missingFirst;
                TopoDS_Vertex missingSecond;
                if (getEdgeVerticesOriented(
                        bestMissingEdge,
                        missingFirst,
                        missingSecond)) {
                    const gp_Pnt openFirstPoint = BRep_Tool::Pnt(openFirst);
                    const gp_Pnt openSecondPoint = BRep_Tool::Pnt(openSecond);
                    const double directDistance = std::max(
                        openFirstPoint.Distance(BRep_Tool::Pnt(missingFirst)),
                        openSecondPoint.Distance(BRep_Tool::Pnt(missingSecond)));
                    const double reversedDistance = std::max(
                        openFirstPoint.Distance(BRep_Tool::Pnt(missingSecond)),
                        openSecondPoint.Distance(BRep_Tool::Pnt(missingFirst)));
                    const bool reverseMissingVertices =
                        reversedDistance < directDistance;
                    const double endpointDistance = std::min(
                        directDistance,
                        reversedDistance);
                    const double endpointMergeTolerance = std::max({
                        Precision::Confusion(),
                        myParams.maxEndpointTolerance,
                        BRep_Tool::Tolerance(openFirst),
                        BRep_Tool::Tolerance(openSecond),
                        BRep_Tool::Tolerance(missingFirst),
                        BRep_Tool::Tolerance(missingSecond) });

                    TopoDS_Edge supplementalDefectiveEdge;
                    double supplementalDefectiveEndpointDistance =
                        std::numeric_limits<double>::max();
                    for (const TopoDS_Edge& removedEdge : removedEdges) {
                        TopoDS_Vertex removedFirst;
                        TopoDS_Vertex removedSecond;
                        if (!getEdgeVerticesOriented(
                                removedEdge,
                                removedFirst,
                                removedSecond)) {
                            continue;
                        }
                        const double directRemovedDistance = std::max(
                            openFirstPoint.Distance(BRep_Tool::Pnt(removedFirst)),
                            openSecondPoint.Distance(BRep_Tool::Pnt(removedSecond)));
                        const double reversedRemovedDistance = std::max(
                            openFirstPoint.Distance(BRep_Tool::Pnt(removedSecond)),
                            openSecondPoint.Distance(BRep_Tool::Pnt(removedFirst)));
                        const double removedEndpointDistance = std::min(
                            directRemovedDistance,
                            reversedRemovedDistance);
                        if (removedEndpointDistance <
                            supplementalDefectiveEndpointDistance) {
                            supplementalDefectiveEndpointDistance =
                                removedEndpointDistance;
                            supplementalDefectiveEdge = removedEdge;
                        }
                    }

                    // 被删除的独占 edge 本来就是错误边界，补充边不需要复现它。
                    // 能否继续使用原 surface，应由补充边到该 surface 的实测偏差决定。
                    const double supplementalDeviationBudget = std::max({
                        Precision::Confusion(),
                        myParams.maxBoundaryToSurfaceDeviation,
                        BRep_Tool::Tolerance(choice.faceToRebuild) });
                    bool supplementalBoundaryIsValid = false;
                    if (endpointDistance <= endpointMergeTolerance &&
                        supplementalDefectiveEndpointDistance <=
                            endpointMergeTolerance &&
                        !supplementalDefectiveEdge.IsNull()) {
                        supplementalBoundaryIsValid =
                            measureEdgeToFaceSupportDeviation(
                                bestMissingEdge,
                                choice.faceToRebuild,
                                std::max(8, candidate.usedSampleCount),
                                supplementalSupportDeviation,
                                supplementalFailureReason) &&
                            supplementalSupportDeviation <=
                                supplementalDeviationBudget;
                        if (!supplementalBoundaryIsValid &&
                            supplementalFailureReason.empty()) {
                            std::ostringstream failure;
                            failure << "support deviation exceeds budget: deviation="
                                    << supplementalSupportDeviation
                                    << ", maximum="
                                    << supplementalDeviationBudget;
                            supplementalFailureReason = failure.str();
                        }
                    }

                    if (supplementalBoundaryIsValid) {
                        BRepTools_ReShape vertexReshaper;
                        const auto registerVertexMerge = [&vertexReshaper](
                            const TopoDS_Vertex& sourceVertex,
                            const TopoDS_Vertex& targetVertex) {
                            if (sourceVertex.IsSame(targetVertex)) {
                                return;
                            }
                            TopoDS_Vertex sourceKey = sourceVertex;
                            sourceKey.Orientation(TopAbs_FORWARD);
                            TopoDS_Vertex replacementValue = targetVertex;
                            replacementValue.Orientation(TopAbs_FORWARD);
                            vertexReshaper.Replace(sourceKey, replacementValue);
                        };
                        registerVertexMerge(
                            missingFirst,
                            reverseMissingVertices ? openSecond : openFirst);
                        registerVertexMerge(
                            missingSecond,
                            reverseMissingVertices ? openFirst : openSecond);

                        Candidate supplementalCandidate;
                        supplementalCandidate.edgeA =
                            supplementalDefectiveEdge;
                        supplementalCandidate.edgeB = bestMissingEdge;
                        supplementalCandidate.faceA = choice.faceToRebuild;
                        supplementalCandidate.faceB = bestMissingEdgeOwner;
                        supplementalCandidate.lengthA =
                            GeometryUtils::ComputeEdgeLength(
                                supplementalDefectiveEdge);
                        supplementalCandidate.lengthB =
                            GeometryUtils::ComputeEdgeLength(bestMissingEdge);
                        const double supplementalMaximumLength = std::max(
                            supplementalCandidate.lengthA,
                            supplementalCandidate.lengthB);
                        if (supplementalMaximumLength > Precision::Confusion()) {
                            supplementalCandidate.lengthRatio = std::min(
                                supplementalCandidate.lengthA,
                                supplementalCandidate.lengthB) /
                                supplementalMaximumLength;
                        }

                        std::vector<Candidate> retryCandidates = trialCandidates;
                        retryCandidates.push_back(supplementalCandidate);
                        const TopoDS_Shape mergedTrialShape =
                            vertexReshaper.Apply(trialShape);
                        const auto mapCandidateShape = [&vertexReshaper](
                            const TopoDS_Shape& sourceShape,
                            TopAbs_ShapeEnum expectedType,
                            TopoDS_Shape& mappedShape) {
                            mappedShape = vertexReshaper.Apply(sourceShape);
                            if (mappedShape.IsNull() ||
                                mappedShape.ShapeType() != expectedType) {
                                return false;
                            }
                            mappedShape.Orientation(sourceShape.Orientation());
                            return true;
                        };
                        bool mappingIsValid = true;
                        for (Candidate& retryCandidate : retryCandidates) {
                            const Candidate sourceCandidate = retryCandidate;
                            TopoDS_Shape mappedShape;
                            mappingIsValid = mapCandidateShape(
                                sourceCandidate.edgeA,
                                TopAbs_EDGE,
                                mappedShape);
                            if (mappingIsValid) {
                                retryCandidate.edgeA = TopoDS::Edge(mappedShape);
                                mappingIsValid = mapCandidateShape(
                                    sourceCandidate.edgeB,
                                    TopAbs_EDGE,
                                    mappedShape);
                            }
                            if (mappingIsValid) {
                                retryCandidate.edgeB = TopoDS::Edge(mappedShape);
                                mappingIsValid = mapCandidateShape(
                                    sourceCandidate.faceA,
                                    TopAbs_FACE,
                                    mappedShape);
                            }
                            if (mappingIsValid) {
                                retryCandidate.faceA = TopoDS::Face(mappedShape);
                                mappingIsValid = mapCandidateShape(
                                    sourceCandidate.faceB,
                                    TopAbs_FACE,
                                    mappedShape);
                            }
                            if (mappingIsValid) {
                                retryCandidate.faceB = TopoDS::Face(mappedShape);
                            }
                            if (!mappingIsValid) {
                                break;
                            }
                        }

                        std::string mergeValidationError;
                        if (!mergedTrialShape.IsNull() && mappingIsValid &&
                            validateWholeShape(
                                trialShape,
                                mergedTrialShape,
                                TopologyUtils::CountFreeEdges(trialShape),
                                mergeValidationError)) {
                            appendMessage(
                                messages_,
                                "Crossed-face merged geometrically coincident boundary vertices "
                                "to expose a supplemental canonical edge");
                            std::string retryError;
                            if (tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
                                    mergedTrialShape,
                                    retryCandidates,
                                    repairedShape,
                                    retryError,
                                    false)) {
                                return true;
                            }
                            errorMessage = "retry after boundary vertex merge failed: "
                                + retryError;
                            return false;
                        }
                    }
                }
            }

            std::ostringstream message;
            message << "surviving boundary is not a degree-two graph"
                    << ", wireCount=" << defectiveWires.size()
                    << ", candidatePairs=" << groupedChoices.size()
                    << ", boundaryEdges=" << unorderedBoundary.size()
                    << ", removedEdges=" << removedEdges.size()
                    << ", invalidVertices=" << invalidVertexCount
                    << ", openVertices=" << openVertexIndices.size()
                    << ", maximumVertexDegree=" << maximumVertexDegree;
            if (std::isfinite(nearestOpenVertexDistance)) {
                message << ", nearestOpenVertexDistance="
                        << nearestOpenVertexDistance;
            }
            if (std::isfinite(bestMissingEdgeEndpointDistance)) {
                message << ", bestMissingEdgeIndex=" << bestMissingEdgeIndex
                        << ", bestMissingEdgeOwnerCount="
                        << bestMissingEdgeOwnerCount
                        << ", bestMissingEdgeEndpointDistance="
                        << bestMissingEdgeEndpointDistance;
            }
            if (supplementalSupportDeviation > 0.0 ||
                !supplementalFailureReason.empty()) {
                message << ", supplementalSupportDeviation="
                        << supplementalSupportDeviation;
            }
            if (!supplementalFailureReason.empty()) {
                message << ", supplementalFailure="
                        << supplementalFailureReason;
            }
            errorMessage = message.str();
            return false;
        }

        // 删除后可能留下外环和内环。每次从一条未使用 edge 开始，按顶点
        // 连通关系构造一个独立闭环，直到所有 edge 都被归入某个 wire。
        struct BoundaryLoop {
            std::vector<TopoDS_Edge> edges;
            std::vector<TopoDS_Face> owners;
        };
        std::vector<bool> usedLinks(links.size(), false);
        std::vector<BoundaryLoop> boundaryLoops;
        std::vector<TopoDS_Edge> boundaryEdges;
        std::vector<TopoDS_Face> boundaryOwners;
        while (true) {
            const auto unusedLink = std::find(
                usedLinks.begin(),
                usedLinks.end(),
                false);
            if (unusedLink == usedLinks.end()) {
                break;
            }

            const std::size_t firstLinkIndex = static_cast<std::size_t>(
                std::distance(usedLinks.begin(), unusedLink));
            const TopoDS_Vertex loopStart = links[firstLinkIndex].firstVertex;
            TopoDS_Vertex currentVertex = loopStart;
            std::size_t nextLinkIndex = firstLinkIndex;
            BoundaryLoop loop;
            for (std::size_t step = 0; step < links.size(); ++step) {
                if (nextLinkIndex >= links.size() || usedLinks[nextLinkIndex]) {
                    errorMessage = "surviving boundary loop reused an edge";
                    return false;
                }

                const BoundaryLink& link = links[nextLinkIndex];
                TopoDS_Edge orientedEdge = link.edge;
                TopoDS_Vertex nextVertex;
                if (link.firstVertex.IsSame(currentVertex)) {
                    nextVertex = link.lastVertex;
                }
                else if (link.lastVertex.IsSame(currentVertex)) {
                    orientedEdge.Reverse();
                    nextVertex = link.firstVertex;
                }
                else {
                    errorMessage = "surviving boundary edge is not incident to the current vertex";
                    return false;
                }
                loop.edges.push_back(orientedEdge);
                loop.owners.push_back(link.ownerFace);
                usedLinks[nextLinkIndex] = true;
                currentVertex = nextVertex;
                if (currentVertex.IsSame(loopStart)) {
                    break;
                }

                const int currentVertexIndex =
                    boundaryVertices.FindIndex(currentVertex);
                if (currentVertexIndex <= 0) {
                    errorMessage = "cannot locate the current surviving boundary vertex";
                    return false;
                }
                nextLinkIndex = links.size();
                for (const std::size_t incidentLinkIndex :
                     incidentLinks[static_cast<std::size_t>(currentVertexIndex)]) {
                    if (!usedLinks[incidentLinkIndex]) {
                        if (nextLinkIndex != links.size()) {
                            errorMessage = "surviving boundary has more than one next edge";
                            return false;
                        }
                        nextLinkIndex = incidentLinkIndex;
                    }
                }
                if (nextLinkIndex == links.size()) {
                    errorMessage = "surviving boundary ended before its loop closed";
                    return false;
                }
            }
            if (!currentVertex.IsSame(loopStart) || loop.edges.size() < 3) {
                errorMessage = "surviving boundary contains an open or undersized loop";
                return false;
            }
            boundaryEdges.insert(
                boundaryEdges.end(),
                loop.edges.begin(),
                loop.edges.end());
            boundaryOwners.insert(
                boundaryOwners.end(),
                loop.owners.begin(),
                loop.owners.end());
            boundaryLoops.push_back(std::move(loop));
        }
        if (boundaryLoops.empty()) {
            errorMessage = "surviving boundary contains no closed loop";
            return false;
        }

        // 删除坏 face 后只重组边界，不拟合新曲面。新 face 继续使用坏 face
        // 原来的 support surface，这与 CADFix 中基于已有 surface 重建一致。
        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) originalSurface =
            BRep_Tool::Surface(forwardDefectiveFace, surfaceLocation);
        if (originalSurface.IsNull()) {
            errorMessage = "defective face has no original support surface";
            return false;
        }

        const double existingFaceTolerance = std::max(
            Precision::Confusion(),
            BRep_Tool::Tolerance(forwardDefectiveFace));
        const double faceTolerance = std::min(
            std::max(existingFaceTolerance, myParams.faceTol),
            myParams.maxFaceTol);
        std::vector<TopoDS_Wire> boundaryWires;
        boundaryWires.reserve(boundaryLoops.size());
        for (const BoundaryLoop& loop : boundaryLoops) {
            BRepBuilderAPI_MakeWire wireMaker;
            for (const TopoDS_Edge& boundaryEdge : loop.edges) {
                wireMaker.Add(boundaryEdge);
            }
            if (!wireMaker.IsDone()) {
                errorMessage = "failed to construct a surviving DeleteAndRebuildOnOriginalSurface wire";
                return false;
            }

            TopoDS_Wire boundaryWire = wireMaker.Wire();
            int sameDirectionCount = 0;
            int oppositeDirectionCount = 0;
            for (std::size_t edgeIndex = 0;
                 edgeIndex < loop.edges.size();
                 ++edgeIndex) {
                TopAbs_Orientation ownerOrientation = TopAbs_FORWARD;
                std::string orientationError;
                if (!findUniqueEdgeUseOrientation(
                        loop.owners[edgeIndex],
                        loop.edges[edgeIndex],
                        ownerOrientation,
                        orientationError)) {
                    errorMessage = "cannot orient surviving boundary wire: "
                        + orientationError;
                    return false;
                }
                if (loop.edges[edgeIndex].Orientation() == ownerOrientation) {
                    ++sameDirectionCount;
                }
                else {
                    ++oppositeDirectionCount;
                }
            }
            if (sameDirectionCount > 0 && oppositeDirectionCount > 0) {
                errorMessage = "surviving faces require inconsistent orientation within one patch wire";
                return false;
            }
            if (sameDirectionCount > 0) {
                boundaryWire.Reverse();
            }
            boundaryWires.push_back(boundaryWire);
        }

        BRep_Builder builder;
        TopoDS_Face patchFace;
        builder.MakeFace(
            patchFace,
            originalSurface,
            surfaceLocation,
            faceTolerance);
        builder.NaturalRestriction(
            patchFace,
            BRep_Tool::NaturalRestriction(forwardDefectiveFace));
        patchFace.Orientation(TopAbs_FORWARD);
        for (const TopoDS_Wire& boundaryWire : boundaryWires) {
            builder.Add(patchFace, boundaryWire);
        }

        double maximumAppliedTolerance = pcurveTolerance;
        double maximumBoundaryDeviation = 0.0;
        for (TopoDS_Edge& boundaryEdge : boundaryEdges) {
            double edgeDeviation = 0.0;
            double edgeTolerance = 0.0;
            if (!projectEdgeToFaceSurface(
                    boundaryEdge,
                    patchFace,
                    myParams.maxBoundaryToSurfaceDeviation,
                    myParams.maxRebuiltEdgeTolerance,
                    edgeDeviation,
                    edgeTolerance,
                    errorMessage)) {
                errorMessage = "failed to prepare DeleteAndRebuildOnOriginalSurface boundary pcurve: "
                    + errorMessage;
                return false;
            }
            maximumBoundaryDeviation = std::max(
                maximumBoundaryDeviation,
                edgeDeviation);
            maximumAppliedTolerance = std::max(
                maximumAppliedTolerance,
                edgeTolerance);
        }
        builder.UpdateFace(
            patchFace,
            faceTolerance);

        // 新 patch 的 face use 方向由所有 surviving face 的共享边方向共同决定。
        int sameDirectionCount = 0;
        int oppositeDirectionCount = 0;
        for (std::size_t edgeIndex = 0;
             edgeIndex < boundaryEdges.size();
             ++edgeIndex) {
            TopAbs_Orientation patchOrientation = TopAbs_FORWARD;
            TopAbs_Orientation ownerOrientation = TopAbs_FORWARD;
            std::string orientationError;
            if (!findUniqueEdgeUseOrientation(
                    patchFace,
                    boundaryEdges[edgeIndex],
                    patchOrientation,
                    orientationError) ||
                !findUniqueEdgeUseOrientation(
                    boundaryOwners[edgeIndex],
                    boundaryEdges[edgeIndex],
                    ownerOrientation,
                    orientationError)) {
                errorMessage = "cannot determine DeleteAndRebuildOnOriginalSurface boundary orientation: "
                    + orientationError;
                return false;
            }
            if (patchOrientation == ownerOrientation) {
                ++sameDirectionCount;
            }
            else {
                ++oppositeDirectionCount;
            }
        }
        if (sameDirectionCount > 0 && oppositeDirectionCount > 0) {
            errorMessage = "surviving faces require inconsistent patch orientations";
            return false;
        }
        if (sameDirectionCount > 0) {
            patchFace.Reverse();
        }

        BRepCheck_Analyzer patchAnalyzer(patchFace, Standard_True);
        if (!patchAnalyzer.IsValid()) {
            bool orientationCorrected = false;
            if (boundaryWires.size() <= 8) {
                const std::size_t combinationCount =
                    static_cast<std::size_t>(1) << boundaryWires.size();
                for (std::size_t mask = 1;
                     mask < combinationCount;
                     ++mask) {
                    BRep_Builder alternativeBuilder;
                    TopoDS_Face alternativeFace;
                    alternativeBuilder.MakeFace(
                        alternativeFace,
                        originalSurface,
                        surfaceLocation,
                        faceTolerance);
                    alternativeBuilder.NaturalRestriction(
                        alternativeFace,
                        BRep_Tool::NaturalRestriction(forwardDefectiveFace));
                    alternativeFace.Orientation(TopAbs_FORWARD);
                    for (std::size_t wireIndex = 0;
                         wireIndex < boundaryWires.size();
                         ++wireIndex) {
                        TopoDS_Wire alternativeWire = boundaryWires[wireIndex];
                        if ((mask & (static_cast<std::size_t>(1) << wireIndex)) != 0) {
                            alternativeWire.Reverse();
                        }
                        alternativeBuilder.Add(alternativeFace, alternativeWire);
                    }
                    alternativeBuilder.UpdateFace(
                        alternativeFace,
                        faceTolerance);

                    BRepCheck_Analyzer alternativeAnalyzer(
                        alternativeFace,
                        Standard_True);
                    if (alternativeAnalyzer.IsValid()) {
                        patchFace = alternativeFace;
                        orientationCorrected = true;
                        std::ostringstream orientationMessage;
                        orientationMessage
                            << "Crossed-face patch wire orientation combination selected: mask="
                            << mask;
                        appendMessage(messages_, orientationMessage.str());
                        break;
                    }
                }
            }

            if (!orientationCorrected) {
                // 组合仍无效时，最后只让 ShapeFix_Face 整理外环/内环方向。
                ShapeFix_Face faceFix(patchFace);
                faceFix.SetPrecision(faceTolerance);
                faceFix.SetMaxTolerance(std::max(
                    faceTolerance,
                    myParams.maxFaceTol));
                faceFix.FixWireMode() = 0;
                faceFix.FixOrientationMode() = 1;
                faceFix.FixMissingSeamMode() = 0;
                faceFix.FixAddNaturalBoundMode() = 0;
                faceFix.Perform();

                const TopoDS_Shape fixedShape = faceFix.Result();
                if (!fixedShape.IsNull() &&
                    fixedShape.ShapeType() == TopAbs_FACE) {
                    patchFace = TopoDS::Face(fixedShape);
                }
                else {
                    patchFace = faceFix.Face();
                }
                if (patchFace.IsNull()) {
                    errorMessage = "ShapeFix_Face returned a null crossed-face patch";
                    return false;
                }

                BRepCheck_Analyzer fixedPatchAnalyzer(patchFace, Standard_True);
                if (!fixedPatchAnalyzer.IsValid()) {
                    errorMessage = "DeleteAndRebuildOnOriginalSurface patch is invalid after local orientation fix: "
                        + RepairDiagnostics::DescribeBRepCheckFailures(
                            fixedPatchAnalyzer,
                            patchFace);
                    return false;
                }
                appendMessage(
                    messages_,
                    "Crossed-face patch wire orientations corrected by local ShapeFix_Face");
            }
        }

        // wire 方向组合或 ShapeFix_Face 可能改变 edge use 方向，因此必须在
        // 最终 patch 上重新协调一次与所有 surviving owner face 的相对方向。
        int finalSameDirectionCount = 0;
        int finalOppositeDirectionCount = 0;
        for (std::size_t edgeIndex = 0;
             edgeIndex < boundaryEdges.size();
             ++edgeIndex) {
            TopAbs_Orientation patchEdgeOrientation = TopAbs_FORWARD;
            TopAbs_Orientation ownerEdgeOrientation = TopAbs_FORWARD;
            std::string orientationError;
            if (!findUniqueEdgeUseOrientation(
                    patchFace,
                    boundaryEdges[edgeIndex],
                    patchEdgeOrientation,
                    orientationError) ||
                !findUniqueEdgeUseOrientation(
                    boundaryOwners[edgeIndex],
                    boundaryEdges[edgeIndex],
                    ownerEdgeOrientation,
                    orientationError)) {
                errorMessage = "cannot coordinate final crossed-face patch orientation: "
                    + orientationError;
                return false;
            }
            if (patchEdgeOrientation == ownerEdgeOrientation) {
                ++finalSameDirectionCount;
            }
            else {
                ++finalOppositeDirectionCount;
            }
        }
        if (finalSameDirectionCount > 0 &&
            finalOppositeDirectionCount > 0) {
            errorMessage = "final crossed-face patch edges require inconsistent face orientations";
            return false;
        }
        if (finalSameDirectionCount > 0) {
            patchFace.Reverse();
        }

        const int freeEdgesBefore = TopologyUtils::CountFreeEdges(trialShape);
        TopoDS_Face sourceKey = choice.faceToRebuild;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Face replacementValue = patchFace;
        const TopAbs_Orientation sourceOrientation = choice.faceToRebuild.Orientation();
        const TopAbs_Orientation patchOrientation = patchFace.Orientation();
        if ((sourceOrientation != TopAbs_FORWARD &&
             sourceOrientation != TopAbs_REVERSED) ||
            (patchOrientation != TopAbs_FORWARD &&
             patchOrientation != TopAbs_REVERSED)) {
            errorMessage = "defective face or DeleteAndRebuild patch has unsupported orientation";
            return false;
        }
        replacementValue.Orientation(
            sourceOrientation == patchOrientation
                ? TopAbs_FORWARD
                : TopAbs_REVERSED);
        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        repairedShape = reshaper.Apply(trialShape);
        if (repairedShape.IsNull() ||
            repairedShape.ShapeType() != trialShape.ShapeType()) {
            errorMessage = "DeleteAndRebuildOnOriginalSurface changed the root shape type";
            repairedShape.Nullify();
            return false;
        }
        // 提交前先保持原有的整体 BRepCheck 增量和 free-edge 约束。
        if (!validateWholeShape(
                trialShape,
                repairedShape,
                freeEdgesBefore,
                errorMessage)) {
            repairedShape.Nullify();
            return false;
        }

        TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
        TopologyUtils::BuildEdgeFaceMap(repairedShape, repairedEdgeFaceMap);
        TopTools_ListOfShape removedEdgeFaces;
        if (findFacesForSameEdge(
                repairedEdgeFaceMap,
                choice.edgeToReplace,
                removedEdgeFaces)) {
            errorMessage = "defective edge is still used after DeleteAndRebuildOnOriginalSurface";
            repairedShape.Nullify();
            return false;
        }
        for (const TopoDS_Edge& removedEdge : removedEdges) {
            TopTools_ListOfShape removedFaces;
            if (findFacesForSameEdge(
                    repairedEdgeFaceMap,
                    removedEdge,
                    removedFaces)) {
                errorMessage = "an edge owned only by the defective face survived DeleteAndRebuildOnOriginalSurface";
                repairedShape.Nullify();
                return false;
            }
        }

        TopTools_IndexedMapOfShape resultFaces;
        TopExp::MapShapes(repairedShape, TopAbs_FACE, resultFaces);
        for (int faceIndex = 1; faceIndex <= resultFaces.Extent(); ++faceIndex) {
            if (resultFaces.FindKey(faceIndex).IsSame(choice.faceToRebuild)) {
                errorMessage = "defective face is still present after DeleteAndRebuildOnOriginalSurface";
                repairedShape.Nullify();
                return false;
            }
        }

        for (const TopoDS_Edge& boundaryEdge : boundaryEdges) {
            TopTools_ListOfShape adjacentFaces;
            if (!findFacesForSameEdge(
                    repairedEdgeFaceMap,
                    boundaryEdge,
                    adjacentFaces) ||
                adjacentFaces.Extent() != 2) {
                errorMessage = "DeleteAndRebuildOnOriginalSurface boundary edge is not shared by two faces";
                repairedShape.Nullify();
                return false;
            }

            TopAbs_Orientation orientations[2] = {
                TopAbs_FORWARD,
                TopAbs_FORWARD
            };
            int faceIndex = 0;
            for (TopTools_ListIteratorOfListOfShape faceIterator(adjacentFaces);
                 faceIterator.More();
                 faceIterator.Next(), ++faceIndex) {
                const TopoDS_Face adjacentFace = TopoDS::Face(faceIterator.Value());
                std::string orientationError;
                if (!findUniqueEdgeUseOrientation(
                        adjacentFace,
                        boundaryEdge,
                        orientations[faceIndex],
                        orientationError) ||
                    !GeometryUtils::HasValidPCurveOnFace(
                        boundaryEdge,
                        adjacentFace,
                        maximumAppliedTolerance,
                        &orientationError)) {
                    errorMessage = "DeleteAndRebuildOnOriginalSurface boundary validation failed: "
                        + orientationError;
                    repairedShape.Nullify();
                    return false;
                }
            }
            if (orientations[0] == orientations[1]) {
                errorMessage = "DeleteAndRebuildOnOriginalSurface boundary edge has the same orientation in both faces";
                repairedShape.Nullify();
                return false;
            }
        }

        std::ostringstream message;
        message << "Crossed-face DeleteAndRebuildOnOriginalSurface accepted: selection="
                << selectionMessage
                << ", missingWire="
                << defectiveWireError
                << ", boundaryEdges=" << boundaryEdges.size()
                << ", removedEdges=" << removedEdges.size()
                << ", maxBoundaryDeviation=" << maximumBoundaryDeviation
                << ", freeEdgesBefore=" << freeEdgesBefore
                << ", freeEdgesAfter=" << TopologyUtils::CountFreeEdges(repairedShape);
        appendMessage(messages_, message.str());
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while deleting and rebuilding crossed face on its original surface: "
            + failureText(failure);
        repairedShape.Nullify();
        return false;
    }
    return true;
}

HealNarrowGap::OriginalSurfaceRepairStatus
HealNarrowGap::tryRepairCrossedFaceOnOriginalSurface(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    double& maxBoundaryDeviation,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    maxBoundaryDeviation = 0.0;

    TopoDS_Shape trialShape;
    Candidate trialCandidate;
    if (!createIsolatedCandidateTrial(
            currentShape,
            candidate,
            trialShape,
            trialCandidate,
            errorMessage)) {
        return OriginalSurfaceRepairStatus::RejectCandidate;
    }

    EdgeChoice choice;
    std::string selectionMessage;
    if (!chooseCrossedFaceBoundaryEdgeChoice(
            trialCandidate,
            choice,
            selectionMessage,
            errorMessage)) {
        return OriginalSurfaceRepairStatus::RejectCandidate;
    }
    TopoDS_Wire correctedWire;
    std::vector<TopoDS_Edge> boundaryEdges;
    if (!buildCorrectedFaceBoundary(
            choice,
            correctedWire,
            boundaryEdges,
            errorMessage)) {
        return OriginalSurfaceRepairStatus::RejectCandidate;
    }

    double fixedPCurveTolerance = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.keptEdge),
        BRep_Tool::Tolerance(choice.fixedFace),
        myParams.pcurveTol });
    std::string fixedPCurveError;
    if (!GeometryUtils::HasValidPCurveOnFace(
            choice.keptEdge,
            choice.fixedFace,
            fixedPCurveTolerance,
            &fixedPCurveError)) {
        errorMessage = "canonical kept edge is not a trusted boundary: invalid pcurve on fixed face: "
            + fixedPCurveError;
        return OriginalSurfaceRepairStatus::RejectCandidate;
    }

    try {
        TopoDS_Face forwardDefectiveFace = choice.faceToRebuild;
        forwardDefectiveFace.Orientation(TopAbs_FORWARD);
        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) surface =
            BRep_Tool::Surface(forwardDefectiveFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "defective face has no original support surface";
            return OriginalSurfaceRepairStatus::RejectCandidate;
        }

        const double existingFaceTolerance = std::max(
            Precision::Confusion(),
            BRep_Tool::Tolerance(forwardDefectiveFace));
        const double faceTolerance = std::min(
            std::max(existingFaceTolerance, myParams.faceTol),
            myParams.maxFaceTol);

        BRep_Builder builder;
        TopoDS_Face rebuiltFace;
        builder.MakeFace(
            rebuiltFace,
            surface,
            surfaceLocation,
            faceTolerance);
        builder.NaturalRestriction(
            rebuiltFace,
            BRep_Tool::NaturalRestriction(forwardDefectiveFace));
        rebuiltFace.Orientation(TopAbs_FORWARD);
        builder.Add(rebuiltFace, correctedWire);

        TopoDS_Edge replacementUse;
        for (TopoDS_Edge& boundaryEdge : boundaryEdges) {
            if (boundaryEdge.IsSame(choice.keptEdge)) {
                replacementUse = boundaryEdge;
                break;
            }
        }
        if (replacementUse.IsNull()) {
            errorMessage = "corrected boundary does not contain the canonical kept edge";
            return OriginalSurfaceRepairStatus::RejectCandidate;
        }

        double appliedTolerance = 0.0;
        if (!projectEdgeToFaceSurface(
                replacementUse,
                rebuiltFace,
                myParams.maxBoundaryToSurfaceDeviation,
                myParams.maxRebuiltEdgeTolerance,
                maxBoundaryDeviation,
                appliedTolerance,
                errorMessage)) {
            return OriginalSurfaceRepairStatus::NeedsSurfaceRebuild;
        }
        builder.UpdateFace(
            rebuiltFace,
            std::max(faceTolerance, appliedTolerance));
        rebuiltFace.Orientation(choice.faceToRebuild.Orientation());

        BRepCheck_Analyzer faceAnalyzer(rebuiltFace, Standard_True);
        if (!faceAnalyzer.IsValid()) {
            errorMessage = "corrected wire is valid but cannot form a valid face on the original surface: "
                + RepairDiagnostics::DescribeBRepCheckFailures(
                    faceAnalyzer,
                    rebuiltFace);
            return OriginalSurfaceRepairStatus::NeedsSurfaceRebuild;
        }

        if (!replaceCrossedFaceAndValidate(
                trialShape,
                choice,
                rebuiltFace,
                std::max(appliedTolerance, fixedPCurveTolerance),
                repairedShape,
                errorMessage)) {
            return OriginalSurfaceRepairStatus::RejectCandidate;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while rebuilding crossed face on original surface: "
            + failureText(failure);
        return OriginalSurfaceRepairStatus::NeedsSurfaceRebuild;
    }

    return OriginalSurfaceRepairStatus::Success;
}

bool HealNarrowGap::tryRepairCrossedFaceWithRebuiltSurface(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& rebuildMethod,
    double& maxBoundaryDeviation,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    rebuildMethod.clear();
    maxBoundaryDeviation = 0.0;

    const auto attemptMethod = [&](
        bool usePlanarSurface,
        TopoDS_Shape& methodResult,
        double& methodDeviation,
        std::string& methodError) {
        TopoDS_Shape trialShape;
        Candidate trialCandidate;
        if (!createIsolatedCandidateTrial(
                currentShape,
                candidate,
                trialShape,
                trialCandidate,
                methodError)) {
            return false;
        }

        EdgeChoice choice;
        std::string selectionMessage;
        if (!chooseCrossedFaceBoundaryEdgeChoice(
                trialCandidate,
                choice,
                selectionMessage,
                methodError)) {
            return false;
        }
        TopoDS_Wire correctedWire;
        std::vector<TopoDS_Edge> boundaryEdges;
        if (!buildCorrectedFaceBoundary(
                choice,
                correctedWire,
                boundaryEdges,
                methodError)) {
            return false;
        }

        double fixedPCurveTolerance = std::max(
            BRep_Tool::Tolerance(choice.keptEdge),
            myParams.pcurveTol);
        std::string fixedPCurveError;
        if (!GeometryUtils::HasValidPCurveOnFace(
                choice.keptEdge,
                choice.fixedFace,
                fixedPCurveTolerance,
                &fixedPCurveError)) {
            methodError = "canonical kept edge is not a trusted boundary: invalid pcurve on fixed face: "
                + fixedPCurveError;
            return false;
        }

        try {
            const double faceTolerance = std::min(
                std::max(
                    BRep_Tool::Tolerance(choice.faceToRebuild),
                    myParams.faceTol),
                myParams.maxFaceTol);
            TopoDS_Face patchFace;
            BRep_Builder builder;

            if (usePlanarSurface) {
                BRepBuilderAPI_FindPlane planeFinder(
                    correctedWire,
                    myParams.maxBoundaryToSurfaceDeviation);
                if (!planeFinder.Found() || planeFinder.Plane().IsNull()) {
                    methodError = "corrected boundary is not planar within the deviation budget";
                    return false;
                }
                BRepBuilderAPI_MakeFace faceMaker(
                    planeFinder.Plane()->Pln(),
                    correctedWire,
                    Standard_True);
                if (!faceMaker.IsDone() || faceMaker.Face().IsNull()) {
                    methodError = "failed to build a planar crossed-face patch";
                    return false;
                }
                patchFace = faceMaker.Face();
            }
            else {
                constexpr int fillingDegree = 3;
                constexpr int fillingPointsOnCurve = 15;
                constexpr int fillingIterations = 2;
                constexpr int fillingMaxDegree = 6;
                constexpr int fillingMaxSegments = 16;
                constexpr int maxPolesPerDirection = 80;
                constexpr int maxControlPoints = 4096;

                BRepFill_Filling filling;
                filling.SetResolParam(
                    fillingDegree,
                    fillingPointsOnCurve,
                    fillingIterations,
                    Standard_False);
                filling.SetConstrParam(
                    myParams.pcurveTol,
                    myParams.maxBoundaryToSurfaceDeviation,
                    1.0e-2,
                    1.0e-1);
                filling.SetApproxParam(
                    fillingMaxDegree,
                    fillingMaxSegments);
                for (const TopoDS_Edge& boundaryEdge : boundaryEdges) {
                    filling.Add(boundaryEdge, GeomAbs_C0, Standard_True);
                }
                filling.Build();
                if (!filling.IsDone() || filling.Face().IsNull()) {
                    methodError = "BRepFill_Filling did not build a crossed-face patch";
                    return false;
                }
                const double g0Error = filling.G0Error();
                if (!std::isfinite(g0Error) ||
                    g0Error > myParams.maxBoundaryToSurfaceDeviation) {
                    std::ostringstream message;
                    message << "Filling G0 error exceeds the boundary budget: error="
                            << g0Error
                            << ", maximum="
                            << myParams.maxBoundaryToSurfaceDeviation;
                    methodError = message.str();
                    return false;
                }

                const Handle(Geom_Surface) fillingSurface =
                    BRep_Tool::Surface(filling.Face());
                const Handle(Geom_BSplineSurface) bsplineSurface =
                    Handle(Geom_BSplineSurface)::DownCast(fillingSurface);
                if (bsplineSurface.IsNull()) {
                    methodError = "Filling did not produce a B-spline support surface";
                    return false;
                }
                if (bsplineSurface->UDegree() > fillingMaxDegree ||
                    bsplineSurface->VDegree() > fillingMaxDegree ||
                    bsplineSurface->NbUPoles() > maxPolesPerDirection ||
                    bsplineSurface->NbVPoles() > maxPolesPerDirection ||
                    bsplineSurface->NbUPoles() * bsplineSurface->NbVPoles()
                        > maxControlPoints) {
                    methodError = "Filling support surface exceeds the fixed complexity budget";
                    return false;
                }

                builder.MakeFace(
                    patchFace,
                    fillingSurface,
                    TopLoc_Location(),
                    faceTolerance);
                patchFace.Orientation(TopAbs_FORWARD);
                builder.Add(patchFace, correctedWire);
            }

            methodDeviation = 0.0;
            double maximumAppliedTolerance = fixedPCurveTolerance;
            for (TopoDS_Edge& boundaryEdge : boundaryEdges) {
                double edgeDeviation = 0.0;
                double edgeTolerance = 0.0;
                if (!projectEdgeToFaceSurface(
                        boundaryEdge,
                        patchFace,
                        myParams.maxBoundaryToSurfaceDeviation,
                        myParams.maxRebuiltEdgeTolerance,
                        edgeDeviation,
                        edgeTolerance,
                        methodError)) {
                    return false;
                }
                methodDeviation = std::max(methodDeviation, edgeDeviation);
                maximumAppliedTolerance = std::max(
                    maximumAppliedTolerance,
                    edgeTolerance);
            }
            builder.UpdateFace(
                patchFace,
                std::max(faceTolerance, maximumAppliedTolerance));
            patchFace.Orientation(choice.faceToRebuild.Orientation());

            BRepCheck_Analyzer faceAnalyzer(patchFace, Standard_True);
            if (!faceAnalyzer.IsValid()) {
                methodError = "rebuilt support surface produced an invalid face: "
                    + RepairDiagnostics::DescribeBRepCheckFailures(
                        faceAnalyzer,
                        patchFace);
                return false;
            }

            if (!replaceCrossedFaceAndValidate(
                    trialShape,
                    choice,
                    patchFace,
                    maximumAppliedTolerance,
                    methodResult,
                    methodError)) {
                return false;
            }
        }
        catch (const Standard_Failure& failure) {
            methodError = "OCCT failure while rebuilding crossed-face support surface: "
                + failureText(failure);
            return false;
        }
        return true;
    };

    std::string planarError;
    if (attemptMethod(
            true,
            repairedShape,
            maxBoundaryDeviation,
            planarError)) {
        rebuildMethod = "Planar";
        return true;
    }

    std::string fillingError;
    if (attemptMethod(
            false,
            repairedShape,
            maxBoundaryDeviation,
            fillingError)) {
        rebuildMethod = "Filling";
        return true;
    }

    errorMessage = "Level 2 surface rebuild failed:"
        "\n  [Planar]\n"
        + RepairDiagnostics::IndentMultiline(
            RepairDiagnostics::MakeFailureReasonReadable(planarError),
            "    ")
        + "\n  [Filling]\n"
        + RepairDiagnostics::IndentMultiline(
            RepairDiagnostics::MakeFailureReasonReadable(fillingError),
            "    ");
    return false;
}

bool HealNarrowGap::tryRepairByRebuildingCrossedFaceBoundary(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();

    EdgeChoice choice;
    std::string selectionMessage;
    if (!chooseCrossedFaceBoundaryEdgeChoice(
            candidate,
            choice,
            selectionMessage,
            errorMessage)) {
        return false;
    }
    appendMessage(
        messages_,
        "Crossed-face rebuild selection: " + selectionMessage);

    // 先确认存在 crossed-face 证据，再判断两条候选边是否足够接近。
    // 面本身已经声明的 tolerance 高于配置值时，至少尊重输入模型的已有精度。
    const double faceBasedGapTolerance = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.faceA),
        BRep_Tool::Tolerance(candidate.faceB) });
    const double originalSurfaceGapLimit = std::max(
        myParams.crossedFaceOriginalSurfaceMaxGap,
        faceBasedGapTolerance);
    if (candidate.maxGap > originalSurfaceGapLimit) {
        std::ostringstream message;
        message << "crossed-face gap requires a fitted support surface, which is not enabled: maxGap="
                << candidate.maxGap
                << ", originalSurfaceGapLimit=" << originalSurfaceGapLimit
                << ", configuredLimit="
                << myParams.crossedFaceOriginalSurfaceMaxGap
                << ", faceBasedTolerance=" << faceBasedGapTolerance;
        errorMessage = message.str();
        return false;
    }

    std::string defectiveWireError;
    if (!hasTraversableCrossedFaceBoundary(
            choice.faceToRebuild,
            choice.edgeToReplace,
            defectiveWireError)) {
        appendMessage(
            messages_,
            "Crossed-face boundary mode: DeleteAndRebuildOnOriginalSurface, reason="
                + defectiveWireError);
        if (tryRepairByDeletingAndRebuildingCrossedFaceOnOriginalSurface(
                currentShape,
                candidate,
                repairedShape,
                errorMessage)) {
            return true;
        }
        errorMessage = "DeleteAndRebuildOnOriginalSurface failed: " + errorMessage;
        return false;
    }

    double level1Deviation = 0.0;
    std::string level1Error;
    const OriginalSurfaceRepairStatus level1Status =
        tryRepairCrossedFaceOnOriginalSurface(
            currentShape,
            candidate,
            repairedShape,
            level1Deviation,
            level1Error);

    if (level1Status == OriginalSurfaceRepairStatus::Success) {
        std::ostringstream message;
        message << "Crossed-face boundary Level 1: status=Success"
                << ", maxBoundaryDeviation=" << level1Deviation;
        appendMessage(messages_, message.str());
        appendMessage(
            messages_,
            "Candidate repair accepted by crossed-face boundary rebuild: level=1, method=OriginalSurface");
        return true;
    }

    if (level1Status == OriginalSurfaceRepairStatus::RejectCandidate) {
        appendMessage(
            messages_,
            "Crossed-face boundary Level 1: status=RejectCandidate");
        errorMessage = "Level 1 rejected candidate: " + level1Error;
        return false;
    }

    appendMessage(
        messages_,
        "Crossed-face boundary Level 1: status=NeedsSurfaceRebuild");
    appendMessage(
        messages_,
        "Crossed-face fitted-surface fallback: enabled=false");
    errorMessage = "original support surface cannot contain the corrected boundary; "
        "fitted-surface rebuild is not enabled: " + level1Error;
    return false;
}

bool HealNarrowGap::buildSurfaceIntersectionEdgeForCandidate(
    const Candidate& candidate,
    double endpointDistanceTolerance,
    TopoDS_Edge& intersectionEdge,
    std::string& errorMessage) const
{
    intersectionEdge.Nullify();
    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "invalid candidate for surface intersection";
        return false;
    }

    try {
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        if (!getEdgeVerticesOriented(candidate.edgeA, aStart, aEnd) ||
            !getEdgeVerticesOriented(candidate.edgeB, bStart, bEnd)) {
            errorMessage = "failed to get oriented edge vertices for surface intersection";
            return false;
        }

        if (!((aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) ||
              (aStart.IsSame(bEnd) && aEnd.IsSame(bStart)))) {
            errorMessage = "candidate edges do not share endpoints for surface intersection";
            return false;
        }

        BRepAlgoAPI_Section section(candidate.faceA, candidate.faceB, Standard_False);
        if (myParams.surfaceIntersectionFuzzyTol > 0.0) {
            section.SetFuzzyValue(myParams.surfaceIntersectionFuzzyTol);
        }
        section.Approximation(Standard_True);
        section.ComputePCurveOn1(Standard_True);
        section.ComputePCurveOn2(Standard_True);
        section.Build();
        if (!section.IsDone() || section.HasErrors()) {
            errorMessage = "BRepAlgoAPI_Section failed";
            return false;
        }

        const TopoDS_Shape sectionShape = section.Shape();
        if (sectionShape.IsNull()) {
            errorMessage = "surface intersection returned null shape";
            return false;
        }

        BRepAdaptor_Curve edgeCurveA(candidate.edgeA);
        BRepAdaptor_Curve edgeCurveB(candidate.edgeB);
        const gp_Pnt startPoint = BRep_Tool::Pnt(aStart);
        const gp_Pnt endPoint = BRep_Tool::Pnt(aEnd);
        const double endpointTol = std::max({
            Precision::Confusion(),
            endpointDistanceTolerance,
            BRep_Tool::Tolerance(aStart),
            BRep_Tool::Tolerance(aEnd)
            });
        const int sampleCount = std::max(8, std::min(32, computeSampleCountForCandidate(candidate)));

        bool found = false;
        double bestScore = std::numeric_limits<double>::max();
        TopoDS_Edge bestEdge;
        int sectionEdgeCount = 0;
        int rejectedByEndpointDistance = 0;
        int rejectedByGapDistance = 0;
        int rejectedByMakeEdge = 0;
        double bestRejectedStartDistance = std::numeric_limits<double>::max();
        double bestRejectedEndDistance = std::numeric_limits<double>::max();
        double bestRejectedEndpointDistance = std::numeric_limits<double>::max();
        double bestRejectedGapDistance = std::numeric_limits<double>::max();

        for (TopExp_Explorer edgeExp(sectionShape, TopAbs_EDGE); edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge sectionEdge = TopoDS::Edge(edgeExp.Current());
            if (sectionEdge.IsNull() || isDegeneratedEdge(sectionEdge) || !has3DCurve(sectionEdge)) {
                continue;
            }

            double first = 0.0;
            double last = 0.0;
            Handle(Geom_Curve) curve = curveForEdgeInWorld(sectionEdge, first, last);
            if (curve.IsNull()) {
                continue;
            }
            ++sectionEdgeCount;

            double startParam = first;
            double endParam = last;
            double startDistance = 0.0;
            double endDistance = 0.0;
            if (!projectPointOnCurve(startPoint, curve, first, last, startParam, startDistance) ||
                !projectPointOnCurve(endPoint, curve, first, last, endParam, endDistance)) {
                continue;
            }

            if (std::abs(endParam - startParam) <= Precision::PConfusion()) {
                continue;
            }

            // 构造 edge 前先校验原端点是否真正落在 section curve 上。
            // 否则 MakeEdge 可能成功，但后续 pcurve 校验会在端点或近端采样处失败。
            const double endpointDistance = std::max(startDistance, endDistance);
            if (endpointDistance > endpointTol) {
                ++rejectedByEndpointDistance;
                if (endpointDistance < bestRejectedEndpointDistance) {
                    bestRejectedEndpointDistance = endpointDistance;
                    bestRejectedStartDistance = startDistance;
                    bestRejectedEndDistance = endDistance;
                }
                continue;
            }

            double maxDistanceToGapEdges = 0.0;
            double meanDistanceToGapEdges = 0.0;
            bool distanceFailed = false;
            for (int sample = 0; sample < sampleCount; ++sample) {
                const double t = static_cast<double>(sample) / static_cast<double>(sampleCount - 1);
                const double param = startParam + (endParam - startParam) * t;
                const gp_Pnt point = curve->Value(param);

                double distanceToA = 0.0;
                double distanceToB = 0.0;
                if (!nearestDistanceToEdge(point, edgeCurveA, distanceToA) ||
                    !nearestDistanceToEdge(point, edgeCurveB, distanceToB)) {
                    distanceFailed = true;
                    break;
                }

                const double distanceToGap = std::min(distanceToA, distanceToB);
                maxDistanceToGapEdges = std::max(maxDistanceToGapEdges, distanceToGap);
                meanDistanceToGapEdges += distanceToGap;
            }

            if (distanceFailed) {
                continue;
            }
            meanDistanceToGapEdges /= static_cast<double>(sampleCount);

            // section edge 必须靠近原 gap 区域；否则复杂 face-face 相交时可能选到其他交线。
            if (maxDistanceToGapEdges > myParams.surfaceIntersectionMaxGap) {
                ++rejectedByGapDistance;
                bestRejectedGapDistance = std::min(bestRejectedGapDistance, maxDistanceToGapEdges);
                continue;
            }

            Handle(Geom_Curve) edgeCurve = curve;
            double edgeStartParam = startParam;
            double edgeEndParam = endParam;
            if (edgeEndParam < edgeStartParam) {
                // 保证新 edge 的底层曲线方向和候选起点/终点一致，后续 face use 再沿用原 wire 方向。
                edgeCurve = curve->Reversed();
                edgeStartParam = curve->ReversedParameter(startParam);
                edgeEndParam = curve->ReversedParameter(endParam);
            }

            BRepBuilderAPI_MakeEdge edgeMaker(edgeCurve, aStart, aEnd, edgeStartParam, edgeEndParam);
            if (!edgeMaker.IsDone()) {
                // 如果原顶点不在交线上，说明交线不能在不重建相邻边的情况下直接用于当前 wire。
                ++rejectedByMakeEdge;
                continue;
            }

            const double endpointScore = std::max(startDistance, endDistance);
            const double score =
                endpointScore * 4.0 +
                meanDistanceToGapEdges +
                maxDistanceToGapEdges * 0.25;
            if (score < bestScore) {
                bestScore = score;
                bestEdge = edgeMaker.Edge();
                found = true;
            }
        }

        if (!found || bestEdge.IsNull()) {
            std::ostringstream os;
            os << "no usable surface intersection edge close to candidate gap"
               << ": sectionEdges=" << sectionEdgeCount
               << ", endpointDistanceTolerance=" << endpointTol
               << ", rejectedByEndpointDistance=" << rejectedByEndpointDistance
               << ", rejectedByGapDistance=" << rejectedByGapDistance
               << ", rejectedByMakeEdge=" << rejectedByMakeEdge;
            if (bestRejectedEndpointDistance < std::numeric_limits<double>::max()) {
                os << ", bestRejectedStartDistance=" << bestRejectedStartDistance
                   << ", bestRejectedEndDistance=" << bestRejectedEndDistance;
            }
            if (bestRejectedGapDistance < std::numeric_limits<double>::max()) {
                os << ", bestRejectedGapDistance=" << bestRejectedGapDistance;
            }
            errorMessage = os.str();
            return false;
        }

        intersectionEdge = bestEdge;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "surface intersection edge construction failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealNarrowGap::tryRepairBySurfaceIntersection(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }
    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "Invalid candidate for surface intersection";
        return false;
    }
    if (candidate.faceA.IsSame(candidate.faceB)) {
        errorMessage = "surface intersection does not support two edges on the same face";
        return false;
    }
    if (candidate.maxGap > myParams.surfaceIntersectionMaxGap) {
        std::ostringstream os;
        os << "surface intersection gap too large: maxGap=" << candidate.maxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        errorMessage = os.str();
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    // SurfaceIntersection 产生的新 edge 在替换后使用原两条 edge 的最大 tolerance 校验 pcurve。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.edgeA),
        BRep_Tool::Tolerance(candidate.edgeB)
        });

    const double existingFaceTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.faceA),
        BRep_Tool::Tolerance(candidate.faceB)
        });
    const double faceTol = std::min(
        std::max(existingFaceTol, myParams.faceTol),
        myParams.maxFaceTol);

    TopoDS_Edge intersectionEdge;
    if (!buildSurfaceIntersectionEdgeForCandidate(
            candidate,
            myParams.maxEndpointTolerance,
            intersectionEdge,
            errorMessage)) {
        return false;
    }

    TopoDS_Edge edgeAInWire;
    if (!findEdgeUseInFaceWire(candidate.faceA, candidate.edgeA, edgeAInWire, errorMessage)) {
        errorMessage = "edgeA not found in faceA wire for surface intersection: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeBInWire;
    if (!findEdgeUseInFaceWire(candidate.faceB, candidate.edgeB, edgeBInWire, errorMessage)) {
        errorMessage = "edgeB not found in faceB wire for surface intersection: " + errorMessage;
        return false;
    }

    // 新 edge 的使用方向跟随被替换 edge 在各自 wire 中的方向，保持 face 内侧判断依赖的 wire 方向关系。
    TopoDS_Edge intersectionUseA = intersectionEdge;
    intersectionUseA.Orientation(edgeAInWire.Orientation());
    TopoDS_Edge intersectionUseB = intersectionEdge;
    intersectionUseB.Orientation(edgeBInWire.Orientation());

    if (!ensurePCurveOnFace(intersectionUseA, candidate.faceA, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add surface intersection pcurve on faceA: " + errorMessage;
        return false;
    }
    if (!ensurePCurveOnFace(intersectionUseB, candidate.faceB, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add surface intersection pcurve on faceB: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceA;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceA,
            candidate.edgeA,
            intersectionUseA,
            faceTol,
            rebuiltFaceA,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceA with surface intersection edge: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceB;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceB,
            candidate.edgeB,
            intersectionUseB,
            faceTol,
            rebuiltFaceB,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceB with surface intersection edge: " + errorMessage;
        return false;
    }

    BRepCheck_Analyzer analyzerFaceA(rebuiltFaceA);
    if (!analyzerFaceA.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceA is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceA, rebuiltFaceA));
    }
    BRepCheck_Analyzer analyzerFaceB(rebuiltFaceB);
    if (!analyzerFaceB.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceB is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceB, rebuiltFaceB));
    }

    TopoDS_Shape trialShape;
    try {
        BRepTools_ReShape reshaper;
        reshaper.Replace(candidate.faceA, rebuiltFaceA);
        reshaper.Replace(candidate.faceB, rebuiltFaceB);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "Replaced shape is null after surface intersection";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Surface intersection shape replacement failed: " + failureText(failure);
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape intersectionFaces;
    if (!findFacesForSameEdge(repairedEdgeFaceMap, intersectionEdge, intersectionFaces)) {
        errorMessage = "surface intersection edge not found after shape replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(intersectionFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge intersectionUseInFace;
        if (!findEdgeUseInFaceWire(face, intersectionEdge, intersectionUseInFace, pcurveRepairError)) {
            errorMessage = "Failed to find surface intersection edge use after replacement: " + pcurveRepairError;
            return false;
        }

        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                intersectionUseInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "Failed to repair surface intersection pcurve after replacement: "
                + pcurveRepairError;
            return false;
        }
    }

    if (myParams.runShapeFixAfterRepair) {
        try {
            const double fixPrecision = std::max(Precision::Confusion(), myParams.shapeFixPrecision);
            const double fixMaxTol = std::max(fixPrecision, myParams.shapeFixMaxTolerance);
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(fixPrecision);
            fixer.SetMaxTolerance(fixMaxTol);
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix failed after surface intersection: " + failureText(failure);
            return false;
        }
    }

    if (!validateEdgeReplacement(
            currentShape,
            trialShape,
            intersectionEdge,
            freeEdgesBefore,
            pcurveTol,
            errorMessage)) {
        errorMessage = "surface intersection validation failed: " + errorMessage;
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by surface intersection: candidateMaxGap=" << myParams.candidateMaxGap
           << ", repairMaxGap=" << myParams.repairMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", surfaceIntersectionFuzzyTol=" << myParams.surfaceIntersectionFuzzyTol
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidate.maxGap=" << candidate.maxGap
           << ", candidate.meanGap=" << candidate.meanGap
           << ", sampleCount=" << candidate.usedSampleCount
           << ", inheritedPCurveTol=" << pcurveTol
           << ", faceTol=" << faceTol
           << ", shapeFixPrecision=" << myParams.shapeFixPrecision
           << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
        appendMessage(messages_, os.str());
    }

    repairedShape = trialShape;
    return true;
}

bool HealNarrowGap::buildMidlineEdgeForCandidate(
    const Candidate& candidate,
    TopoDS_Edge& midlineEdge,
    std::string& errorMessage) const
{
    midlineEdge.Nullify();
    if (candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "candidate edge is null";
        return false;
    }

    try {
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        if (!getEdgeVerticesOriented(candidate.edgeA, aStart, aEnd) ||
            !getEdgeVerticesOriented(candidate.edgeB, bStart, bEnd)) {
            errorMessage = "failed to get oriented edge vertices for midline";
            return false;
        }

        TopoDS_Vertex bParamStartVertex;
        TopoDS_Vertex bParamEndVertex;
        if (aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) {
            bParamStartVertex = bStart;
            bParamEndVertex = bEnd;
        }
        else if (aStart.IsSame(bEnd) && aEnd.IsSame(bStart)) {
            bParamStartVertex = bEnd;
            bParamEndVertex = bStart;
        }
        else {
            errorMessage = "candidate edges do not share endpoints for midline";
            return false;
        }

        const double uA0 = BRep_Tool::Parameter(aStart, candidate.edgeA);
        const double uA1 = BRep_Tool::Parameter(aEnd, candidate.edgeA);
        const double uB0 = BRep_Tool::Parameter(bParamStartVertex, candidate.edgeB);
        const double uB1 = BRep_Tool::Parameter(bParamEndVertex, candidate.edgeB);

        BRepAdaptor_Curve curveA(candidate.edgeA);
        BRepAdaptor_Curve curveB(candidate.edgeB);
        const int sampleCount = computeSampleCountForCandidate(candidate);

        std::vector<gp_Pnt> midPoints;
        midPoints.reserve(static_cast<std::size_t>(sampleCount));
        for (int i = 0; i < sampleCount; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(sampleCount - 1);
            gp_Pnt point;
            if (i == 0) {
                // 端点沿用原拓扑顶点，保证新 edge 能与原 wire 中相邻 edge 精确连接。
                point = BRep_Tool::Pnt(aStart);
            }
            else if (i + 1 == sampleCount) {
                point = BRep_Tool::Pnt(aEnd);
            }
            else {
                const gp_Pnt pointA = curveA.Value(uA0 + (uA1 - uA0) * t);
                const gp_Pnt pointB = curveB.Value(uB0 + (uB1 - uB0) * t);
                point.SetCoord(
                    0.5 * (pointA.X() + pointB.X()),
                    0.5 * (pointA.Y() + pointB.Y()),
                    0.5 * (pointA.Z() + pointB.Z()));
            }

            if (!midPoints.empty() &&
                point.Distance(midPoints.back()) <= Precision::Confusion()) {
                if (i + 1 == sampleCount) {
                    errorMessage = "midline edge has zero-length end segment";
                    return false;
                }
                continue;
            }
            midPoints.push_back(point);
        }

        if (midPoints.size() < 2) {
            errorMessage = "midline edge has insufficient points";
            return false;
        }

        if (midPoints.size() == 2) {
            BRepBuilderAPI_MakeEdge edgeMaker(aStart, aEnd);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create straight midline edge";
                return false;
            }
            midlineEdge = edgeMaker.Edge();
        }
        else {
            TColgp_Array1OfPnt points(1, static_cast<Standard_Integer>(midPoints.size()));
            for (std::size_t i = 0; i < midPoints.size(); ++i) {
                points.SetValue(static_cast<Standard_Integer>(i + 1), midPoints[i]);
            }

            GeomAPI_PointsToBSpline splineBuilder(points);
            Handle(Geom_BSplineCurve) midlineCurve = splineBuilder.Curve();
            if (midlineCurve.IsNull()) {
                errorMessage = "failed to interpolate midline curve";
                return false;
            }

            BRepBuilderAPI_MakeEdge edgeMaker(midlineCurve, aStart, aEnd);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create midline edge";
                return false;
            }
            midlineEdge = edgeMaker.Edge();
        }

        if (midlineEdge.IsNull()) {
            errorMessage = "midline edge is null";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "midline edge construction failed: " + failureText(failure);
        return false;
    }

    return true;
}

bool HealNarrowGap::tryRepairByMidlineCollapse(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }
    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "Invalid candidate for midline collapse";
        return false;
    }
    if (candidate.faceA.IsSame(candidate.faceB)) {
        errorMessage = "midline collapse does not support two edges on the same face";
        return false;
    }
    if (candidate.maxGap > myParams.midlineMaxGap) {
        std::ostringstream os;
        os << "midline gap too large: maxGap=" << candidate.maxGap
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        errorMessage = os.str();
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    // MidlineCollapse 产生的新 edge 在替换后使用原两条 edge 的最大 tolerance 校验 pcurve。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.edgeA),
        BRep_Tool::Tolerance(candidate.edgeB)
        });

    const double existingFaceTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.faceA),
        BRep_Tool::Tolerance(candidate.faceB)
        });
    const double faceTol = std::min(
        std::max(existingFaceTol, myParams.faceTol),
        myParams.maxFaceTol);

    TopoDS_Edge midlineEdge;
    if (!buildMidlineEdgeForCandidate(candidate, midlineEdge, errorMessage)) {
        return false;
    }

    TopoDS_Edge edgeAInWire;
    if (!findEdgeUseInFaceWire(candidate.faceA, candidate.edgeA, edgeAInWire, errorMessage)) {
        errorMessage = "edgeA not found in faceA wire for midline collapse: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeBInWire;
    if (!findEdgeUseInFaceWire(candidate.faceB, candidate.edgeB, edgeBInWire, errorMessage)) {
        errorMessage = "edgeB not found in faceB wire for midline collapse: " + errorMessage;
        return false;
    }

    // 新 edge 的使用方向跟随被替换 edge 在各自 wire 中的方向，而不是重新用几何端点推断。
    TopoDS_Edge midlineUseA = midlineEdge;
    midlineUseA.Orientation(edgeAInWire.Orientation());
    TopoDS_Edge midlineUseB = midlineEdge;
    midlineUseB.Orientation(edgeBInWire.Orientation());

    if (!ensurePCurveOnFace(midlineUseA, candidate.faceA, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add midline pcurve on faceA: " + errorMessage;
        return false;
    }
    if (!ensurePCurveOnFace(midlineUseB, candidate.faceB, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add midline pcurve on faceB: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceA;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceA,
            candidate.edgeA,
            midlineUseA,
            faceTol,
            rebuiltFaceA,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceA with midline edge: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceB;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceB,
            candidate.edgeB,
            midlineUseB,
            faceTol,
            rebuiltFaceB,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceB with midline edge: " + errorMessage;
        return false;
    }

    BRepCheck_Analyzer analyzerFaceA(rebuiltFaceA);
    if (!analyzerFaceA.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceA is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceA, rebuiltFaceA));
    }
    BRepCheck_Analyzer analyzerFaceB(rebuiltFaceB);
    if (!analyzerFaceB.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceB is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceB, rebuiltFaceB));
    }

    TopoDS_Shape trialShape;
    try {
        BRepTools_ReShape reshaper;
        reshaper.Replace(candidate.faceA, rebuiltFaceA);
        reshaper.Replace(candidate.faceB, rebuiltFaceB);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "Replaced shape is null after midline collapse";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Midline shape replacement failed: " + failureText(failure);
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape midlineFaces;
    if (!findFacesForSameEdge(repairedEdgeFaceMap, midlineEdge, midlineFaces)) {
        errorMessage = "midline edge not found after shape replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(midlineFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge midlineUseInFace;
        if (!findEdgeUseInFaceWire(face, midlineEdge, midlineUseInFace, pcurveRepairError)) {
            errorMessage = "Failed to find midline edge use after replacement: " + pcurveRepairError;
            return false;
        }

        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                midlineUseInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "Failed to repair midline pcurve after replacement: " + pcurveRepairError;
            return false;
        }
    }

    if (myParams.runShapeFixAfterRepair) {
        try {
            const double fixPrecision = std::max(Precision::Confusion(), myParams.shapeFixPrecision);
            const double fixMaxTol = std::max(fixPrecision, myParams.shapeFixMaxTolerance);
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(fixPrecision);
            fixer.SetMaxTolerance(fixMaxTol);
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix failed after midline collapse: " + failureText(failure);
            return false;
        }
    }

    if (!validateEdgeReplacement(
            currentShape,
            trialShape,
            midlineEdge,
            freeEdgesBefore,
            pcurveTol,
            errorMessage)) {
        errorMessage = "midline collapse validation failed: " + errorMessage;
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by midline collapse: candidateMaxGap=" << myParams.candidateMaxGap
           << ", repairMaxGap=" << myParams.repairMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidate.maxGap=" << candidate.maxGap
           << ", candidate.meanGap=" << candidate.meanGap
           << ", sampleCount=" << candidate.usedSampleCount
           << ", inheritedPCurveTol=" << pcurveTol
           << ", faceTol=" << faceTol
           << ", shapeFixPrecision=" << myParams.shapeFixPrecision
           << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
        appendMessage(messages_, os.str());
    }

    repairedShape = trialShape;
    return true;
}


bool HealNarrowGap::findEdgeUseInFaceWire(
    const TopoDS_Face& oldFace,
    const TopoDS_Edge& oldEdge,
    TopoDS_Edge& oldEdgeInWire,
    std::string& errorMessage) const
{
    oldEdgeInWire.Nullify();
    if (oldFace.IsNull() || oldEdge.IsNull()) {
        errorMessage = "Invalid face or edge for edge replacement check";
        return false;
    }

    try {
        for (TopExp_Explorer wireExp(oldFace, TopAbs_WIRE); wireExp.More(); wireExp.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExp.Current());
            for (BRepTools_WireExplorer edgeExp(oldWire, oldFace); edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Edge edgeInWire = edgeExp.Current();
                if (!edgeInWire.IsSame(oldEdge)) {
                    continue;
                }

                oldEdgeInWire = edgeInWire;
                return true;
            }
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Failed to find replacement edge in face: " + failureText(failure);
        return false;
    }

    errorMessage = "Old edge not found in face";
    return false;
}

bool HealNarrowGap::rebuildFaceByReplacingEdge(
    const TopoDS_Face& faceToRebuild,
    const TopoDS_Edge& edgeToReplace,
    const TopoDS_Edge& replacementEdgeUse,
    double faceTol,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage) const
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull() || edgeToReplace.IsNull() || replacementEdgeUse.IsNull()) {
        errorMessage = "Invalid face or edge for face to rebuild";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "Replaced face has no support surface";
            return false;
        }

        BRep_Builder builder;
        TopoDS_Face rebuilt;
        builder.MakeFace(rebuilt, surface, surfaceLocation, faceTol);
        builder.NaturalRestriction(rebuilt, BRep_Tool::NaturalRestriction(forwardFace));
        rebuilt.Orientation(TopAbs_FORWARD);

        bool replaced = false;
        for (TopExp_Explorer wireExp(forwardFace, TopAbs_WIRE); wireExp.More(); wireExp.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExp.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            bool wireContainsEdge = false;
            for (BRepTools_WireExplorer edgeExp(oldWire, forwardFace); edgeExp.More(); edgeExp.Next()) {
                const TopoDS_Edge edgeInWire = edgeExp.Current();
                if (edgeInWire.IsSame(edgeToReplace)) {
                    wireContainsEdge = true;
                    break;
                }
            }

            TopoDS_Wire wireToAdd = oldWire;
            if (wireContainsEdge) {
                BRepBuilderAPI_MakeWire wireMaker;
                bool replacedInWire = false;
                for (BRepTools_WireExplorer edgeExp(oldWire, forwardFace); edgeExp.More(); edgeExp.Next()) {
                    const TopoDS_Edge edgeInWire = edgeExp.Current();
                    if (edgeInWire.IsSame(edgeToReplace)) {
                        wireMaker.Add(replacementEdgeUse);
                        replacedInWire = true;
                    }
                    else {
                        wireMaker.Add(edgeInWire);
                    }
                }

                if (!replacedInWire || !wireMaker.IsDone()) {
                    errorMessage = "Failed to rebuild target wire";
                    return false;
                }

                wireToAdd = wireMaker.Wire();
                wireToAdd.Orientation(oldWire.Orientation());

                ShapeFix_Wire wireFix(wireToAdd, forwardFace, faceTol);
                wireFix.ModifyTopologyMode() = Standard_False;
                wireFix.ModifyGeometryMode() = Standard_True;
                wireFix.ClosedWireMode() = Standard_True;
                wireFix.FixReorderMode() = 1;
                wireFix.FixConnectedMode() = 1;
                wireFix.FixEdgeCurvesMode() = 1;
                wireFix.FixReversed2dMode() = 1;
                wireFix.FixSameParameterMode() = 1;
                wireFix.Perform();

                TopoDS_Wire fixedWire = wireFix.WireAPIMake();
                if (fixedWire.IsNull()) {
                    fixedWire = wireFix.Wire();
                }
                if (fixedWire.IsNull()) {
                    errorMessage = "ShapeFix_Wire returned null wire";
                    return false;
                }

                fixedWire.Orientation(oldWire.Orientation());
                wireToAdd = fixedWire;
                replaced = true;
            }

            builder.Add(rebuilt, wireToAdd);
        }

        if (!replaced) {
            errorMessage = "Edge to replace not found while rebuilding face";
            return false;
        }

        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTol);
        faceFix.SetMaxTolerance(std::max(faceTol, myParams.maxFaceTol));
        faceFix.FixWireMode() = 1;
        faceFix.FixOrientationMode() = 1;
        faceFix.FixMissingSeamMode() = 1;
        faceFix.FixAddNaturalBoundMode() = 0;
        faceFix.Perform();

        TopoDS_Shape fixedShape = faceFix.Result();
        if (!fixedShape.IsNull() && fixedShape.ShapeType() == TopAbs_FACE) {
            rebuiltFace = TopoDS::Face(fixedShape);
        }
        else {
            rebuiltFace = faceFix.Face();
        }

        if (rebuiltFace.IsNull()) {
            errorMessage = "ShapeFix_Face returned null face";
            return false;
        }

        rebuiltFace.Orientation(faceToRebuild.Orientation());
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Face rebuild failed: " + failureText(failure);
        return false;
    }

    return true;
}

// 如果 edge 在 face 上缺少 pcurve，则尝试用 ShapeFix_Edge 投影并修正 same-parameter。
bool HealNarrowGap::ensurePCurveOnFace(
    TopoDS_Edge& edge,
    const TopoDS_Face& face,
    double tol,
    std::string& errorMessage) const
{
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "Edge or face is null";
        return false;
    }

    if (!GeometryUtils::RepairPCurveForEdgeOnFace(
            edge,
            face,
            std::max(tol, Precision::Confusion()),
            &errorMessage)) {
        return false;
    }

    return true;
}

// 校验整体 shape 有效，并检查自由边数量没有增加。
bool HealNarrowGap::validateWholeShape(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& repairedShape,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (beforeShape.IsNull()) {
        errorMessage = "Input shape is null";
        return false;
    }

    if (repairedShape.IsNull()) {
        errorMessage = "Repaired shape is null";
        return false;
    }

    std::string existingError;
    const RepairDiagnostics::FailureList existingFailures =
        RepairDiagnostics::CollectBRepCheckFailures(beforeShape, &existingError);
    if (!existingError.empty()) {
        errorMessage = "Existing BRepCheck failure collection failed before replacement: " + existingError;
        return false;
    }

    std::string repairedError;
    const RepairDiagnostics::FailureList repairedFailures =
        RepairDiagnostics::CollectBRepCheckFailures(repairedShape, &repairedError);
    if (!repairedError.empty()) {
        errorMessage = "Repaired shape BRepCheck failed: " + repairedError;
        return false;
    }

    // face/shell 重建后父级 TShape 会变化，不能按对象身份把原有错误误判为新增。
    // 按子形状类型和 status 的数量比较，只拒绝实际增加的错误。
    const RepairDiagnostics::FailureList newFailures =
        RepairDiagnostics::FindBRepCheckFailureCountIncrease(
            existingFailures,
            repairedFailures);

    if (!newFailures.empty()) {
        errorMessage = "New BRepCheck failures after replacement:";
        errorMessage += RepairDiagnostics::DescribeBRepCheckFailures(newFailures);
        if (!existingFailures.empty()) {
            errorMessage += " existing failures:";
            errorMessage += RepairDiagnostics::DescribeBRepCheckFailures(existingFailures);
        }
        return false;
    }

    // freeEdgesAfter 不能大于修复前，否则说明局部替换破坏了边界共享关系。
    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(repairedShape);
    if (freeEdgesAfter > freeEdgesBefore) {
        errorMessage = "free edge count increased";
        return false;
    }

    return true;
}

// 校验 edge 级替换结果：整体合法、自由边未增加、keptEdge 已被多个 face 共享且 pcurve 可用。
bool HealNarrowGap::validateEdgeReplacement(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& repairedShape,
    const TopoDS_Edge& keptEdge,
    int freeEdgesBefore,
    double pcurveValidationTol,
    std::string& errorMessage) const
{
    if (!validateWholeShape(beforeShape, repairedShape, freeEdgesBefore, errorMessage)) {
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(repairedShape, edgeFaceMap);
    TopTools_ListOfShape adjacentFaces;
    if (!findFacesForSameEdge(edgeFaceMap, keptEdge, adjacentFaces)) {
        errorMessage = "Kept edge not found after edge replacement";
        return false;
    }

    if (adjacentFaces.Extent() != 2) {
        errorMessage = "Kept edge is not shared after edge replacement";
        return false;
    }

    for (TopTools_ListIteratorOfListOfShape it(adjacentFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        // 最终 pcurve 校验直接使用替换前两条候选 edge 继承的 tolerance。
        // 这里不再引入 face tolerance、倍率放大或 section fuzzy tolerance。
        const double pcurveCheckTol =
            std::max(pcurveValidationTol, Precision::Confusion());
        std::string localError;
        if (!GeometryUtils::HasValidPCurveOnFace(keptEdge, face, pcurveCheckTol, &localError)) {
            errorMessage = "Kept edge pcurve invalid after edge replacement: " + localError;
            return false;
        }
    }

    return true;
}

// 统一追加日志，便于后续调整日志策略。
void HealNarrowGap::appendMessage(std::vector<std::string>& messages, const std::string& message) const
{
    messages.push_back(message);
}

// Manual verification:
// Positive case: two distinct one-face edges sharing the same vertex pair and
// whose bidirectional nearest gap is below Params::repairMaxGap should yield Perform()>0,
// isModified()==true, two faces sharing the kept TopoDS_Edge, and a valid BRep.
// Negative case: if the bidirectional nearest gap is greater than Params::candidateMaxGap,
// Perform should return 0 and leave isModified()==false.
