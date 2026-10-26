#include "HealClosedNarrowGap.h"

#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

namespace {
constexpr int kMaxAutomaticRepairIterations = 300;
std::string formatIndexList(const std::vector<std::size_t>& indices)
{
    std::ostringstream os;
    os << "[";
    for (std::size_t index = 0; index < indices.size(); ++index) {
        if (index > 0) {
            os << ",";
        }
        os << indices[index];
    }
    os << "]";
    return os.str();
}

} // namespace

HealClosedNarrowGap::HealClosedNarrowGap(const Params& params)
    : params_(params)
{
    params_.maxGapWidth = std::max(params_.maxGapWidth, Precision::Confusion());
    params_.minLengthRatio = GeometryUtils::Clamp01(params_.minLengthRatio);
    params_.minLengthWidthRatio = std::max(1.0, params_.minLengthWidthRatio);
    if (!std::isfinite(params_.maxEndEdgeLengthFactor) ||
        params_.maxEndEdgeLengthFactor <= Precision::Confusion()) {
        params_.maxEndEdgeLengthFactor = 2.0;
    }
    if (!std::isfinite(params_.minAdjacentEdgeLengthRatio) ||
        params_.minAdjacentEdgeLengthRatio < 1.0) {
        params_.minAdjacentEdgeLengthRatio = 3.0;
    }
    params_.widthSampleCount = std::max(2, params_.widthSampleCount);
    params_.minWidthSampleCount = std::max(2, params_.minWidthSampleCount);
    params_.maxWidthSampleCount = std::max(
        params_.minWidthSampleCount,
        params_.maxWidthSampleCount);
    if (!std::isfinite(params_.targetWidthSampleSpacing) ||
        params_.targetWidthSampleSpacing <= Precision::Confusion()) {
        params_.useAdaptiveWidthSampleCount = false;
    }
    params_.sampleCount = std::max(2, params_.sampleCount);
    params_.minSampleCount = std::max(2, params_.minSampleCount);
    params_.maxSampleCount = std::max(params_.minSampleCount, params_.maxSampleCount);
    if (params_.targetSampleSpacing <= Precision::Confusion()) {
        params_.useAdaptiveSampleCount = false;
    }
    params_.midlineRepairTolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
    params_.projectionVertexSnapTolerance = std::max(
        params_.projectionVertexSnapTolerance,
        Precision::Confusion());
    if (!std::isfinite(params_.deleteAndFillTolerance)) {
        params_.deleteAndFillTolerance = 1.0e-3;
    }
    params_.deleteAndFillTolerance = std::max(
        params_.deleteAndFillTolerance,
        Precision::Confusion());
    if (!std::isfinite(params_.deleteAndFillMaxG0Error)) {
        params_.deleteAndFillMaxG0Error = 1.0e-1;
    }
    params_.deleteAndFillMaxG0Error = std::max(
        params_.deleteAndFillMaxG0Error,
        params_.deleteAndFillTolerance);
    params_.deleteAndFillDegree = std::max(2, params_.deleteAndFillDegree);
    params_.deleteAndFillPointsOnCurve = std::max(
        2,
        params_.deleteAndFillPointsOnCurve);
    params_.deleteAndFillIterations = std::max(
        1,
        params_.deleteAndFillIterations);
    params_.deleteAndFillMaxDegree = std::max(
        params_.deleteAndFillDegree,
        params_.deleteAndFillMaxDegree);
    params_.deleteAndFillMaxSegments = std::max(
        1,
        params_.deleteAndFillMaxSegments);
}

void HealClosedNarrowGap::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = shape;
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
}

int HealClosedNarrowGap::perform()
{
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    TopoDS_Shape currentShape = inputShape_;
    const int iterationLimit = computeAutomaticIterationLimit(currentShape);
    {
        std::ostringstream os;
        os << "Automatic repair iteration limit = " << iterationLimit;
        appendMessage(os.str());
    }
    {
        std::ostringstream os;
        os << "Parameters: maxGapWidth=" << params_.maxGapWidth
           << ", minLengthRatio=" << params_.minLengthRatio
           << ", minLengthWidthRatio=" << params_.minLengthWidthRatio
           << ", maxTipAngle=" << params_.maxTipAngle
           << ", detectShortEdgeEnds="
           << (params_.detectShortEdgeEnds ? "true" : "false")
           << ", maxEndEdgeLengthFactor="
           << params_.maxEndEdgeLengthFactor
           << ", minAdjacentEdgeLengthRatio="
           << params_.minAdjacentEdgeLengthRatio
           << ", widthSampleCount=" << params_.widthSampleCount
           << ", useAdaptiveWidthSampleCount="
           << (params_.useAdaptiveWidthSampleCount ? "true" : "false")
           << ", minWidthSampleCount=" << params_.minWidthSampleCount
           << ", maxWidthSampleCount=" << params_.maxWidthSampleCount
           << ", targetWidthSampleSpacing="
           << params_.targetWidthSampleSpacing
           << ", sampleCount=" << params_.sampleCount
           << ", useAdaptiveSampleCount="
           << (params_.useAdaptiveSampleCount ? "true" : "false")
           << ", minSampleCount=" << params_.minSampleCount
           << ", maxSampleCount=" << params_.maxSampleCount
           << ", targetSampleSpacing=" << params_.targetSampleSpacing
           << ", useMidlineCollapse=" << (params_.useMidlineCollapse ? "true" : "false")
           << ", midlineConstructionMode="
           << midlineConstructionModeName(params_.midlineConstructionMode)
           << ", projectionVertexSnapTolerance="
           << params_.projectionVertexSnapTolerance
           << ", useDeleteAndFillDefectiveFace="
           << (params_.useDeleteAndFillDefectiveFace ? "true" : "false")
           << ", deleteAndFillTolerance="
           << params_.deleteAndFillTolerance
           << ", deleteAndFillMaxG0Error="
           << params_.deleteAndFillMaxG0Error
           << ", deleteAndFillDegree="
           << params_.deleteAndFillDegree
           << ", deleteAndFillPointsOnCurve="
           << params_.deleteAndFillPointsOnCurve
           << ", deleteAndFillIterations="
           << params_.deleteAndFillIterations
           << ", deleteAndFillMaxDegree="
           << params_.deleteAndFillMaxDegree
           << ", deleteAndFillMaxSegments="
           << params_.deleteAndFillMaxSegments
           << ", fallbackToSideReplacement=" << (params_.fallbackToSideReplacement ? "true" : "false")
           << ", midlineRepairTolerance=" << params_.midlineRepairTolerance
           << ", runShapeFixAfterRepair="
           << (params_.runShapeFixAfterRepair ? "true" : "false");
        appendMessage(os.str());
    }

    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Original shape", currentShape);

    for (int iter = 0; iter < iterationLimit; ++iter) {
        RepairDiagnostics::AppendIterationHeader(messages_, iter + 1, iter > 0);

        std::vector<ClosedGapCandidate> candidates = extractCandidates(currentShape);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const ClosedGapCandidate& lhs, const ClosedGapCandidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "Candidate count = " << candidates.size();
            appendMessage(os.str());
        }

        if (candidates.empty()) {
            break;
        }

        // 开始尝试本轮候选前只收集一次；失败候选不更新基线，提交后下一轮重建。
        std::string baselineError;
        const RepairDiagnostics::FailureList existingFailures =
            RepairDiagnostics::CollectBRepCheckFailures(currentShape, &baselineError);
        if (!baselineError.empty()) {
            appendMessage("Repair stopped: existing BRepCheck failure collection failed: "
                + baselineError);
            break;
        }

        bool repairedThisIteration = false;
        for (const ClosedGapCandidate& candidate : candidates) {
            TopoDS_Shape repairedShape;
            std::string errorMessage;
            if (!tryRepairCandidate(currentShape, existingFailures, candidate, repairedShape, errorMessage)) {
                std::ostringstream os;
                os << "Candidate failed:"
                   << "\n  chainAEdges=" << candidate.chainA.edges.size()
                   << "\n  chainBEdges=" << candidate.chainB.edges.size()
                   << "\n  startEnd="
                   << gapEndTypeName(candidate.startEnd.type)
                   << "\n  finishEnd="
                   << gapEndTypeName(candidate.finishEnd.type)
                   << "\n  maxWidth=" << candidate.maxWidth
                   << "\n  maxWidthAToB=" << candidate.maxWidthAToB
                   << "\n  maxWidthBToA=" << candidate.maxWidthBToA
                   << "\n  meanWidth=" << candidate.meanWidth
                   << "\n  widthSampleCount="
                   << candidate.usedWidthSampleCount
                   << "\n  lengthRatio=" << candidate.lengthRatio
                   << "\n  lengthWidthRatio=" << candidate.lengthWidthRatio
                   << "\n  endpointSpan=" << candidate.endpointSpan
                   << "\n  endpointSpanWidthRatio="
                   << candidate.endpointSpanWidthRatio
                   << "\n  endpointScore=" << candidate.endpointScore
                   << "\n  confidence=" << candidate.confidence
                   << "\n  reason:"
                   << "\n" << RepairDiagnostics::IndentMultiline(
                       RepairDiagnostics::MakeFailureReasonReadable(errorMessage),
                       "    ");
                appendMessage(os.str());
                continue;
            }

            currentShape = repairedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++repairedCount_;
            repairedThisIteration = true;

            std::ostringstream os;
            os << "Repaired candidate: chainAEdges=" << candidate.chainA.edges.size()
               << ", chainBEdges=" << candidate.chainB.edges.size()
               << ", startEnd="
               << gapEndTypeName(candidate.startEnd.type)
               << ", finishEnd="
               << gapEndTypeName(candidate.finishEnd.type)
               << ", maxWidth=" << candidate.maxWidth
               << ", maxWidthAToB=" << candidate.maxWidthAToB
               << ", maxWidthBToA=" << candidate.maxWidthBToA
               << ", meanWidth=" << candidate.meanWidth
               << ", widthSampleCount="
               << candidate.usedWidthSampleCount
               << ", lengthRatio=" << candidate.lengthRatio
               << ", lengthWidthRatio=" << candidate.lengthWidthRatio
               << ", endpointSpan=" << candidate.endpointSpan
               << ", endpointSpanWidthRatio="
               << candidate.endpointSpanWidthRatio
               << ", endpointScore=" << candidate.endpointScore
               << ", confidence=" << candidate.confidence;
            appendMessage(os.str());
            break;
        }

        if (!repairedThisIteration) {
            appendMessage("Repair stopped: no candidate was submitted in this iteration");
            break;
        }

        if (iter + 1 == iterationLimit) {
            appendMessage("Repair stopped: automatic iteration limit reached");
        }
    }

    resultShape_ = currentShape;
    appendMessage("");
    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Result shape", resultShape_);
    return repairedCount_;
}

int HealClosedNarrowGap::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

bool HealClosedNarrowGap::tryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    if (params_.useMidlineCollapse) {
        std::string midlineError;
        if (tryRepairByMidlineCollapse(currentShape, existingFailures, candidate, repairedShape, midlineError)) {
            return true;
        }

        std::string deleteAndFillError;
        if (params_.useDeleteAndFillDefectiveFace) {
            if (tryRepairByDeletingAndFillingDefectiveFace(
                    currentShape,
                    existingFailures,
                    candidate,
                    repairedShape,
                    deleteAndFillError)) {
                appendMessage(
                    "Warning: repaired by DeleteAndFillDefectiveFace after earlier strategy failed:\n"
                    + RepairDiagnostics::IndentMultiline(
                        RepairDiagnostics::StrategyFailureBlock(
                            "MidlineCollapse",
                            midlineError),
                        "  "));
                return true;
            }
        }
        else {
            deleteAndFillError = "strategy is disabled";
        }

        std::string sideError;
        if (params_.fallbackToSideReplacement) {
            const bool hasShortEdgeEnd =
                candidate.startEnd.type == GapEndType::ShortEdge ||
                candidate.finishEnd.type == GapEndType::ShortEdge;
            const bool asymmetricSourceChains =
                candidate.chainA.edges.size() != candidate.chainB.edges.size();
            if (hasShortEdgeEnd) {
                sideError =
                    "candidate contains a short edge endpoint; "
                    "side replacement cannot collapse and remove endpoint connector edges";
            }
            else if (asymmetricSourceChains) {
                sideError =
                    "candidate uses asymmetric source chains; "
                    "midline collapse is disabled or failed; "
                    "legacy side replacement only supports one-to-one source edge replacement";
            }
            else if (tryRepairBySideReplacement(
                         currentShape,
                         existingFailures,
                         candidate,
                         repairedShape,
                         sideError)) {
                appendMessage(
                    "Warning: repaired by side replacement after earlier strategies failed:\n"
                    + RepairDiagnostics::IndentMultiline(
                        RepairDiagnostics::StrategyFailureBlock(
                            "MidlineCollapse",
                            midlineError),
                        "  ")
                    + "\n"
                    + RepairDiagnostics::IndentMultiline(
                        RepairDiagnostics::StrategyFailureBlock(
                            "DeleteAndFillDefectiveFace",
                            deleteAndFillError),
                        "  "));
                return true;
            }
        }
        else {
            sideError = "strategy is disabled";
        }

        errorMessage = "Closed narrow gap repair failed:\n"
            + RepairDiagnostics::IndentMultiline(
                RepairDiagnostics::StrategyFailureBlock("MidlineCollapse", midlineError), "  ")
            + "\n"
            + RepairDiagnostics::IndentMultiline(
                RepairDiagnostics::StrategyFailureBlock(
                    "DeleteAndFillDefectiveFace",
                    deleteAndFillError),
                "  ")
            + "\n"
            + RepairDiagnostics::IndentMultiline(
                RepairDiagnostics::StrategyFailureBlock("SideReplacement", sideError), "  ");
        return false;
    }

    if (candidate.startEnd.type == GapEndType::ShortEdge ||
        candidate.finishEnd.type == GapEndType::ShortEdge) {
        errorMessage =
            "candidate contains a short edge endpoint; "
            "midline collapse is required to remove endpoint connector edges";
        return false;
    }
    return tryRepairBySideReplacement(
        currentShape,
        existingFailures,
        candidate,
        repairedShape,
        errorMessage);
}

bool HealClosedNarrowGap::tryRepairByMidlineCollapse(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    MidlineConstructionMode resolvedMode =
        MidlineConstructionMode::Auto;
    std::string modeError;
    if (!resolveMidlineConstructionMode(
            candidate,
            resolvedMode,
            modeError)) {
        appendMessage(
            "Midline construction mode: "
            + std::string(midlineConstructionModeName(
                params_.midlineConstructionMode)));
        errorMessage = modeError;
        return false;
    }
    appendMessage(
        "Midline construction mode: "
        + std::string(midlineConstructionModeName(resolvedMode)));

    MidlineChainData midlineData;
    if (!buildMidlineChain(candidate, midlineData, errorMessage)) {
        return false;
    }

    if (midlineData.constructionMode ==
        MidlineConstructionMode::ProjectionCorrespondence) {
        std::ostringstream summary;
        summary << "Projection correspondence midline topology: chainAEdges="
                << candidate.chainA.edges.size()
                << ", chainBEdges=" << candidate.chainB.edges.size()
                << ", chainAVertices=" << candidate.chainA.vertices.size()
                << ", chainBVertices=" << candidate.chainB.vertices.size()
                << ", correspondenceNodeCount="
                << midlineData.correspondenceNodes.size()
                << ", shortEdgeNodeMerges="
                << midlineData.shortEdgeMergeCount
                << ", targetVertices=" << midlineData.chain.vertices.size()
                << ", targetEdges=" << midlineData.chain.edges.size();
        appendMessage(summary.str());

        for (std::size_t nodeIndex = 0;
             nodeIndex < midlineData.correspondenceNodes.size();
             ++nodeIndex) {
            const ProjectionCorrespondenceNode& node =
                midlineData.correspondenceNodes[nodeIndex];
            std::ostringstream nodeLog;
            nodeLog << "correspondence node " << nodeIndex
                    << ": stationA=" << node.stationA
                    << ", stationB=" << node.stationB
                    << ", endpoint="
                    << (node.endpoint ? "true" : "false")
                    << ", sourceVertexA=";
            if (node.hasSourceVertexA) {
                nodeLog << node.sourceVertexIndexA;
            }
            else {
                nodeLog << "none";
            }
            nodeLog << ", sourceVertexB=";
            if (node.hasSourceVertexB) {
                nodeLog << node.sourceVertexIndexB;
            }
            else {
                nodeLog << "none";
            }
            nodeLog << ", snappedPair="
                    << (node.snappedPair ? "true" : "false")
                    << ", mergedForShortEdge="
                    << (node.mergedForShortEdge ? "true" : "false")
                    << ", projectionDistance="
                    << node.projectionDistance;
            appendMessage(nodeLog.str());
        }

        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.replacementsForChainA.size();
             ++edgeIndex) {
            std::ostringstream mapping;
            mapping << "chainA edge " << edgeIndex
                    << " -> target edges "
                    << formatIndexList(
                        midlineData.replacementsForChainA[edgeIndex]);
            appendMessage(mapping.str());
        }
        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.replacementsForChainB.size();
             ++edgeIndex) {
            std::ostringstream mapping;
            mapping << "chainB edge " << edgeIndex
                    << " -> target edges "
                    << formatIndexList(
                        midlineData.replacementsForChainB[edgeIndex]);
            appendMessage(mapping.str());
        }
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    if (!rebuildFacesByCollapsingToMidline(
            currentShape,
            candidate,
            midlineData,
            trialShape,
            errorMessage)) {
        return false;
    }

    appendMessage(
        "Midline adjacency after local rebuild: "
        + RepairDiagnostics::DescribeEdgeFaceCounts(trialShape, midlineData.chain.edges));

    if (params_.runShapeFixAfterRepair) {
        try {
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(params_.midlineRepairTolerance);
            fixer.SetMaxTolerance(std::max(params_.midlineRepairTolerance, params_.maxGapWidth));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();

            const TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix_Shape failed after midline collapse: " + RepairDiagnostics::FailureText(failure);
            return false;
        }
    }

    appendMessage(
        std::string(params_.runShapeFixAfterRepair
            ? "Midline adjacency after final ShapeFix: "
            : "Midline adjacency before validation (final ShapeFix disabled): ")
        + RepairDiagnostics::DescribeEdgeFaceCounts(trialShape, midlineData.chain.edges));

    std::vector<TopoDS_Edge> removedEndpointEdges;
    if (candidate.startEnd.type == GapEndType::ShortEdge) {
        removedEndpointEdges.push_back(candidate.startEnd.connectorEdge);
    }
    if (candidate.finishEnd.type == GapEndType::ShortEdge) {
        removedEndpointEdges.push_back(candidate.finishEnd.connectorEdge);
    }
    if (!validateRepairResult(
            existingFailures,
            trialShape,
            midlineData.chain,
            removedEndpointEdges,
            freeEdgesBefore,
            errorMessage)) {
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by midline collapse: midlineEdges="
           << midlineData.chain.edges.size()
           << ", midlineLength=" << midlineData.chain.length
           << ", repairTolerance=" << params_.midlineRepairTolerance;
        appendMessage(os.str());
    }
    repairedShape = trialShape;
    return true;
}

bool HealClosedNarrowGap::tryRepairBySideReplacement(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    ChainChoice choice = chooseChainToKeep(candidate);
    if (choice.keptChain.edges.empty() || choice.replaceChain.edges.empty()) {
        errorMessage = "invalid chain choice";
        return false;
    }

    if (choice.keptChain.edges.size() != choice.replaceChain.edges.size()) {
        errorMessage =
            "candidate uses asymmetric source chains; "
            "midline collapse is disabled or failed; "
            "legacy side replacement only supports one-to-one source edge replacement";
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    if (!rebuildFacesByReplacingChain(
            currentShape,
            choice.keptChain,
            choice.replaceChain,
            trialShape,
            errorMessage)) {
        return false;
    }

    if (params_.runShapeFixAfterRepair) {
        try {
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(std::max(Precision::Confusion(), params_.maxGapWidth));
            fixer.SetMaxTolerance(std::max(Precision::Confusion(), params_.maxGapWidth * 10.0));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();

            const TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix_Shape failed after closed gap repair: " + RepairDiagnostics::FailureText(failure);
            return false;
        }
    }

    if (!validateRepairResult(
            existingFailures,
            trialShape,
            choice.keptChain,
            std::vector<TopoDS_Edge>(),
            freeEdgesBefore,
            errorMessage)) {
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by side replacement: keptEdges="
           << choice.keptChain.edges.size()
           << ", replacedEdges=" << choice.replaceChain.edges.size()
           << ", keptLength=" << choice.keptChain.length
           << ", replacedLength=" << choice.replaceChain.length;
        appendMessage(os.str());
    }
    repairedShape = trialShape;
    return true;
}

HealClosedNarrowGap::ChainChoice HealClosedNarrowGap::chooseChainToKeep(
    const ClosedGapCandidate& candidate) const
{
    ChainChoice choice;
    bool keepA = true;

    if (params_.preferLongerChainToKeep) {
        keepA = candidate.chainA.length >= candidate.chainB.length;
    }
    else if (candidate.chainA.edges.size() != candidate.chainB.edges.size()) {
        keepA = candidate.chainA.edges.size() <= candidate.chainB.edges.size();
    }

    if (keepA) {
        choice.keptChain = candidate.chainA;
        choice.replaceChain = candidate.chainB;
    }
    else {
        choice.keptChain = candidate.chainB;
        choice.replaceChain = candidate.chainA;
    }
    return choice;
}

bool HealClosedNarrowGap::validateRepairResult(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& afterShape,
    const EdgeChain& keptChain,
    const std::vector<TopoDS_Edge>& removedEdges,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (beforeShape.IsNull() || afterShape.IsNull()) {
        errorMessage = "before or after shape is null";
        return false;
    }

    std::string existingError;
    const RepairDiagnostics::FailureList existingFailures =
        RepairDiagnostics::CollectBRepCheckFailures(beforeShape, &existingError);
    if (!existingError.empty()) {
        errorMessage = "existing BRepCheck failure collection failed: " + existingError;
        return false;
    }

    return validateRepairResult(
        existingFailures, afterShape, keptChain, removedEdges, freeEdgesBefore, errorMessage);
}

// 仅检查本次结果；缓存限于同一轮 currentShape，不缓存任何修复结果。
bool HealClosedNarrowGap::validateRepairResult(
    const RepairDiagnostics::FailureList& existingFailures,
    const TopoDS_Shape& afterShape,
    const EdgeChain& keptChain,
    const std::vector<TopoDS_Edge>& removedEdges,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (afterShape.IsNull()) {
        errorMessage = "before or after shape is null";
        return false;
    }

    std::string repairedError;
    const RepairDiagnostics::FailureList repairedFailures =
        RepairDiagnostics::CollectBRepCheckFailures(afterShape, &repairedError);
    if (!repairedError.empty()) {
        errorMessage = "repaired BRepCheck failure collection failed: " + repairedError;
        return false;
    }

    // face 重建会改变 TShape 身份；按错误类型和状态的数量比较，避免把旧错误迁移误判为新增。
    const RepairDiagnostics::FailureList newFailures =
        RepairDiagnostics::FindBRepCheckFailureCountIncrease(
            existingFailures,
            repairedFailures);

    if (!newFailures.empty() && params_.rejectInvalidResult) {
        errorMessage = "new BRepCheck failures after closed narrow gap repair:"
            + RepairDiagnostics::DescribeBRepCheckFailures(newFailures);
        return false;
    }

    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (freeEdgesAfter > freeEdgesBefore && params_.rejectFreeEdgeIncrease) {
        std::ostringstream os;
        os << "free edge count increased after closed narrow gap repair: before="
           << freeEdgesBefore << ", after=" << freeEdgesAfter;
        errorMessage = os.str();
        return false;
    }

    if (!removedEdges.empty()) {
        TopTools_IndexedMapOfShape resultEdges;
        TopExp::MapShapes(afterShape, TopAbs_EDGE, resultEdges);
        for (const TopoDS_Edge& removedEdge : removedEdges) {
            for (int edgeIndex = 1;
                edgeIndex <= resultEdges.Extent();
                 ++edgeIndex) {
                if (resultEdges.FindKey(edgeIndex).IsSame(removedEdge)) {
                    errorMessage = "an edge scheduled for removal is still present after repair";
                    return false;
                }
            }
        }
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(afterShape, edgeFaceMap);
    for (const TopoDS_Edge& keptEdge : keptChain.edges) {
        int faceCount = 0;
        for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
            const TopoDS_Edge mappedEdge = TopoDS::Edge(edgeFaceMap.FindKey(i));
            if (!mappedEdge.IsNull() && mappedEdge.IsSame(keptEdge)) {
                faceCount += edgeFaceMap.FindFromIndex(i).Extent();
            }
        }

        if (faceCount < 2) {
            errorMessage = "kept chain edge is not shared after repair";
            return false;
        }
    }

    return true;
}

int HealClosedNarrowGap::computeAutomaticIterationLimit(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return 1;
    }

    TopTools_IndexedMapOfShape edgeMap;
    TopExp::MapShapes(shape, TopAbs_EDGE, edgeMap);
    if (edgeMap.Extent() <= 0) {
        return 1;
    }

    return std::min(edgeMap.Extent(), kMaxAutomaticRepairIterations);
}

int HealClosedNarrowGap::computeWidthSampleCount(
    double lengthA,
    double lengthB) const
{
    const int fallbackCount = std::max(2, params_.widthSampleCount);
    if (!params_.useAdaptiveWidthSampleCount ||
        !std::isfinite(params_.targetWidthSampleSpacing) ||
        params_.targetWidthSampleSpacing <= Precision::Confusion()) {
        return fallbackCount;
    }

    if (!std::isfinite(lengthA) ||
        !std::isfinite(lengthB) ||
        lengthA <= Precision::Confusion() ||
        lengthB <= Precision::Confusion()) {
        return fallbackCount;
    }

    const double pairLength = std::max(lengthA, lengthB);
    const double countByLength =
        std::ceil(pairLength / params_.targetWidthSampleSpacing) + 1.0;
    if (!std::isfinite(countByLength) ||
        countByLength >= static_cast<double>(params_.maxWidthSampleCount)) {
        return params_.maxWidthSampleCount;
    }
    if (countByLength <= static_cast<double>(params_.minWidthSampleCount)) {
        return params_.minWidthSampleCount;
    }

    return static_cast<int>(countByLength);
}

int HealClosedNarrowGap::computeMidlineSampleCount(double lengthA, double lengthB) const
{
    const int fallbackCount = std::max(2, params_.sampleCount);
    if (!params_.useAdaptiveSampleCount ||
        params_.targetSampleSpacing <= Precision::Confusion()) {
        return fallbackCount;
    }

    const double pairLength = std::max(lengthA, lengthB);
    if (!std::isfinite(pairLength) || pairLength <= Precision::Confusion()) {
        return fallbackCount;
    }

    // 加一是为了让采样序列同时包含 source edge 的两个端点。
    const int countByLength = static_cast<int>(std::ceil(
        pairLength / params_.targetSampleSpacing)) + 1;
    return std::max(
        params_.minSampleCount,
        std::min(params_.maxSampleCount, countByLength));
}

int HealClosedNarrowGap::findSameEdgeIndex(
    const std::vector<TopoDS_Edge>& edges,
    const TopoDS_Edge& edge) const
{
    if (edge.IsNull()) {
        return -1;
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!edges[i].IsNull() && edges[i].IsSame(edge)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool HealClosedNarrowGap::ensurePCurveOnFace(
    TopoDS_Edge& edge,
    const TopoDS_Face& face,
    double tolerance,
    std::string& errorMessage) const
{
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "edge or face is null";
        return false;
    }

    if (!GeometryUtils::RepairPCurveForEdgeOnFace(
            edge,
            face,
            std::max(tolerance, Precision::Confusion()),
            &errorMessage)) {
        return false;
    }
    return true;
}

void HealClosedNarrowGap::appendMessage(const std::string& message)
{
    messages_.push_back(message);
}
