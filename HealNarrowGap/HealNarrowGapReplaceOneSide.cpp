#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <sstream>

#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_ReShape.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <TopTools_ListOfShape.hxx>

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

bool HealNarrowGap::tryRepairByReplacingOneSide(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
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

    if (candidate.maxGap > myParams.replaceOneSideMaxGap) {
        std::ostringstream os;
        os << "repair gap too large: maxGap=" << candidate.maxGap
           << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
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

    // pcurve 修复和校验至少使用配置精度，并保留两条原 edge 已有的 tolerance。
    // 该容差不参与候选 gap 判断、section 求交或 face/ShapeFix 容差计算。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.keptEdge),
        BRep_Tool::Tolerance(choice.edgeToReplace),
        myParams.pcurveTol
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
        errorMessage = "Shape replacement failed: " + RepairDiagnostics::FailureText(failure);
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
            errorMessage = "ShapeFix failed: " + RepairDiagnostics::FailureText(failure);
            return false;
        }

    }

    if (!validateEdgeReplacement(
            existingFailures,
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
           << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", surfaceIntersectionFuzzyTol=" << myParams.surfaceIntersectionFuzzyTol
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidate.maxGap=" << candidate.maxGap
           << ", candidate.meanGap=" << candidate.meanGap
           << ", sampleCount=" << candidate.usedSampleCount
           << ", effectivePCurveTol=" << pcurveTol
           << ", faceTol=" << faceTol
           << ", shapeFixPrecision=" << myParams.shapeFixPrecision
           << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
        appendMessage(messages_, os.str());
    }

    repairedShape = trialShape;
    return true;
}

