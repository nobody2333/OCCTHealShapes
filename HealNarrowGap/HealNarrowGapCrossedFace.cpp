#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_Filling.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

#include <GeomAbs_Shape.hxx>

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
            + RepairDiagnostics::FailureText(failure);
        return false;
    }
}

// 根据策略尝试提交单个候选修复；不改变主循环“每轮只提交一个候选”的行为。
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
            + RepairDiagnostics::FailureText(failure);
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
                + RepairDiagnostics::FailureText(failure);
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

    // 重建前再次确认两条候选边足够接近，不要求额外的相交证据。
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

