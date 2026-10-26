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
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_ReShape.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <Precision.hxx>
#include <ShapeFix_Face.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <TopTools_ListOfShape.hxx>

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
            + RepairDiagnostics::FailureText(failure);
        repairedShape.Nullify();
        return false;
    }
    return true;
}

