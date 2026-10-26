#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTools_ReShape.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Geometry.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

// 返回 edge 在有向 face 中累计 face orientation 后的唯一 edge-use 方向。
bool HealNarrowGap::findUniqueEdgeUseOrientation(
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

bool HealNarrowGap::containsSameEdge(
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

bool HealNarrowGap::findFacesForSameEdge(const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
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
            + RepairDiagnostics::FailureText(failure);
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
            + RepairDiagnostics::FailureText(failure);
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
            + RepairDiagnostics::FailureText(failure);
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
            + RepairDiagnostics::FailureText(failure);
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
        errorMessage = "Failed to find replacement edge in face: " + RepairDiagnostics::FailureText(failure);
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
        errorMessage = "Face rebuild failed: " + RepairDiagnostics::FailureText(failure);
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

    return validateWholeShape(existingFailures, repairedShape, freeEdgesBefore, errorMessage);
}

// 复用调用方提供的错误基线，不缓存修复结果的 BRepCheck 检查。
bool HealNarrowGap::validateWholeShape(
    const RepairDiagnostics::FailureList& existingFailures,
    const TopoDS_Shape& repairedShape,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (repairedShape.IsNull()) {
        errorMessage = "Repaired shape is null";
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

    return validateSharedEdgeAfterReplacement(
        repairedShape, keptEdge, pcurveValidationTol, errorMessage);
}

// 普通候选使用本轮基线，整体检查标准与中间 trial 的入口保持一致。
bool HealNarrowGap::validateEdgeReplacement(
    const RepairDiagnostics::FailureList& existingFailures,
    const TopoDS_Shape& repairedShape,
    const TopoDS_Edge& keptEdge,
    int freeEdgesBefore,
    double pcurveValidationTol,
    std::string& errorMessage) const
{
    if (!validateWholeShape(existingFailures, repairedShape, freeEdgesBefore, errorMessage)) {
        return false;
    }

    return validateSharedEdgeAfterReplacement(
        repairedShape, keptEdge, pcurveValidationTol, errorMessage);
}

// 整体检查通过后，沿用原有的 adjacency 和 pcurve 校验。
bool HealNarrowGap::validateSharedEdgeAfterReplacement(
    const TopoDS_Shape& repairedShape,
    const TopoDS_Edge& keptEdge,
    double pcurveValidationTol,
    std::string& errorMessage) const
{
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

        // 最终 pcurve 校验使用普通策略的有效容差（原 edge tolerance 与
        // 配置 pcurveTol 的最大值），不引入 face tolerance 或 section fuzzy tolerance。
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
