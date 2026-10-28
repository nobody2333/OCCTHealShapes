#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include <Precision.hxx>
#include <TopoDS.hxx>


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

} // namespace

const char* HealNarrowGap::repairStrategyName(HealNarrowGap::RepairStrategy strategy)
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

// 构造时规范化参数，避免采样数和迭代数出现无效值。
HealNarrowGap::HealNarrowGap(const Params& params)
    : myParams(params)
{
    myParams.candidateMaxGap = std::max(myParams.candidateMaxGap, Precision::Confusion());
    if (myParams.replaceOneSideMaxGap <= 0.0) {
        myParams.replaceOneSideMaxGap = myParams.candidateMaxGap;
    }
    else {
        myParams.replaceOneSideMaxGap =
            std::max(myParams.replaceOneSideMaxGap, Precision::Confusion());
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
    myParams.minLengthRatio = GeometryUtils::Clamp01(myParams.minLengthRatio);
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
        // 本轮候选共用 currentShape 的错误基线；提交后下一轮重新收集。
        std::string baselineError;
        const RepairDiagnostics::FailureList existingFailures =
            RepairDiagnostics::CollectBRepCheckFailures(currentShape, &baselineError);
        if (!baselineError.empty()) {
            appendMessage(messages_,
                "Repair stopped: Existing BRepCheck failure collection failed before replacement: "
                + baselineError);
            break;
        }

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
               << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
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
                    existingFailures,
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
            // 候选，以及进入 Auto 两种常规策略但都失败的候选。
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
                warning << "Info: repaired by CrossedFaceBoundary only after all regular "
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
bool HealNarrowGap::tryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    errorMessage.clear();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }

    const bool sideGapAllowed = candidate.maxGap <= myParams.replaceOneSideMaxGap;
    const bool surfaceIntersectionGapAllowed =
        candidate.maxGap <= myParams.surfaceIntersectionMaxGap;
    const bool midlineGapAllowed = candidate.maxGap <= myParams.midlineMaxGap;

    if (myParams.repairStrategy == RepairStrategy::ReplaceOneSide) {
        return tryRepairByReplacingOneSide(currentShape, existingFailures, candidate, repairedShape, errorMessage);
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

        return tryRepairBySurfaceIntersection(currentShape, existingFailures, candidate, repairedShape, errorMessage);
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

        return tryRepairByMidlineCollapse(currentShape, existingFailures, candidate, repairedShape, errorMessage);
    }

    // Auto 在这里只运行原有常规策略。CrossedFaceBoundary 由 perform() 在
    // 本轮所有常规候选都无法提交后，从普通跳过列表执行后置筛选。
    // 每次失败立即记录，避免后续候选成功时丢失原因；全部失败时仍返回完整汇总。
    std::string autoFailures;
    const auto recordStrategyFailure = [this, &autoFailures](
        const char* strategyName,
        const std::string& reason) {
        const std::string failureBlock =
            RepairDiagnostics::StrategyFailureBlock(strategyName, reason);
        appendMessage(messages_, "Repair strategy failed:\n"
            + RepairDiagnostics::IndentMultiline(failureBlock, "  "));
        if (!autoFailures.empty()) {
            autoFailures += '\n';
        }
        autoFailures += failureBlock;
    };

    std::string sideError;
    if (sideGapAllowed) {
        if (tryRepairByReplacingOneSide(
                currentShape,
                existingFailures,
                candidate,
                repairedShape,
                sideError)) {
            appendMessage(messages_, "Info: repaired by one-side replacement");
            return true;
        }
        else {
            recordStrategyFailure("ReplaceOneSide", sideError);
        }
    }
    else {
        std::ostringstream os;
        os << "repair gap too large for one-side replacement: maxGap=" << candidate.maxGap
           << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        sideError = os.str();
        recordStrategyFailure("ReplaceOneSide", sideError);
    }

    std::string midlineError;
    if (midlineGapAllowed) {
        if (tryRepairByMidlineCollapse(
                currentShape,
                existingFailures,
                candidate,
                repairedShape,
                midlineError)) {
            std::ostringstream info;
            info << "Info: repaired by midline collapse after one-side replacement failed:"
                 << "\n" << RepairDiagnostics::IndentMultiline(autoFailures, "  ");
            appendMessage(messages_, info.str());
            return true;
        }
        else {
            recordStrategyFailure("MidlineCollapse", midlineError);
        }
    }
    else {
        std::ostringstream message;
        message << "midline gap too large: maxGap=" << candidate.maxGap
                << ", midlineMaxGap=" << myParams.midlineMaxGap;
        midlineError = message.str();
        recordStrategyFailure("MidlineCollapse", midlineError);
    }

    errorMessage = "Auto repair failed:\n"
        + RepairDiagnostics::IndentMultiline(autoFailures, "  ");
    return false;
}

// 将一条日志追加到消息列表。
void HealNarrowGap::appendMessage(std::vector<std::string>& messages, const std::string& message) const
{
    messages.push_back(message);
}

// Manual verification:
// Positive case: two distinct one-face edges sharing the same vertex pair and
// whose bidirectional nearest gap is below Params::replaceOneSideMaxGap should yield Perform()>0,
// isModified()==true, two faces sharing the kept TopoDS_Edge, and a valid BRep.
// Negative case: if the bidirectional nearest gap is greater than Params::candidateMaxGap,
// Perform should return 0 and leave isModified()==false.
