#include "HealClosedNarrowGap.h"

#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GeomAbs_Shape.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

namespace {
bool isSameEdgeInVector(const std::vector<TopoDS_Edge>& edges, const TopoDS_Edge& edge)
{
    for (const TopoDS_Edge& existingEdge : edges) {
        if (!existingEdge.IsNull() && existingEdge.IsSame(edge)) {
            return true;
        }
    }
    return false;
}

// 读取 edge 在指定有向 face 中的唯一 edge use 方向。
// TopExp_Explorer 会累计 face orientation，因此该方向可直接用于 shell 邻接关系检查。
bool findUniqueEdgeUseOrientation(
    const TopoDS_Face& face,
    const TopoDS_Edge& edge,
    TopAbs_Orientation& orientation,
    std::string& errorMessage)
{
    if (face.IsNull() || edge.IsNull()) {
        errorMessage = "face or edge is null while checking boundary orientation";
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
            errorMessage =
                "boundary edge use has INTERNAL or EXTERNAL orientation";
            return false;
        }
        orientation = edgeUse.Orientation();
        ++useCount;
    }

    if (useCount != 1) {
        std::ostringstream message;
        message << "boundary edge does not have exactly one use in face: uses="
                << useCount;
        errorMessage = message.str();
        return false;
    }
    return true;
}

} // namespace

bool HealClosedNarrowGap::buildDeleteAndFillPlan(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    DeleteAndFillPlan& plan,
    std::string& errorMessage) const
{
    plan = DeleteAndFillPlan();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    std::vector<TopoDS_Face> relatedFaces;
    auto addRelatedFace = [&](const TopoDS_Face& face) {
        if (!face.IsNull() && !TopologyUtils::IsSameFaceInVector(relatedFaces, face)) {
            relatedFaces.push_back(face);
        }
    };
    for (const TopoDS_Face& face : candidate.chainA.ownerFaces) {
        addRelatedFace(face);
    }
    for (const TopoDS_Face& face : candidate.chainB.ownerFaces) {
        addRelatedFace(face);
    }
    addRelatedFace(candidate.startEnd.connectorOwnerFace);
    addRelatedFace(candidate.finishEnd.connectorOwnerFace);

    auto chainIsEntirelyOwnedBy = [](const EdgeChain& chain, const TopoDS_Face& face) {
        if (chain.edges.empty() || chain.ownerFaces.size() != chain.edges.size()) {
            return false;
        }
        for (const TopoDS_Face& ownerFace : chain.ownerFaces) {
            if (ownerFace.IsNull() || !ownerFace.IsSame(face)) {
                return false;
            }
        }
        return true;
    };
    auto chainUsesFace = [](const EdgeChain& chain, const TopoDS_Face& face) {
        for (const TopoDS_Face& ownerFace : chain.ownerFaces) {
            if (!ownerFace.IsNull() && ownerFace.IsSame(face)) {
                return true;
            }
        }
        return false;
    };

    // 只选择候选相关且确实存在 WireExplorer 遗漏的 face，避免重建正常 face。
    int bestTraversalMismatch = 0;
    bool sourceIsChainA = false;
    for (const TopoDS_Face& relatedFace : relatedFaces) {
        const bool ownsChainA = chainIsEntirelyOwnedBy(candidate.chainA, relatedFace);
        const bool ownsChainB = chainIsEntirelyOwnedBy(candidate.chainB, relatedFace);
        const bool canReplaceChainA = ownsChainA && !chainUsesFace(candidate.chainB, relatedFace);
        const bool canReplaceChainB = ownsChainB && !chainUsesFace(candidate.chainA, relatedFace);
        if (canReplaceChainA == canReplaceChainB) {
            continue;
        }

        TopoDS_Face forwardFace = relatedFace;
        forwardFace.Orientation(TopAbs_FORWARD);
        int traversalMismatch = 0;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire wire = TopoDS::Wire(wireExplorer.Current());
            std::vector<TopoDS_Edge> directUses;
            for (TopoDS_Iterator directIterator(wire, Standard_False, Standard_False);
                 directIterator.More();
                 directIterator.Next()) {
                if (directIterator.Value().ShapeType() == TopAbs_EDGE) {
                    directUses.push_back(TopoDS::Edge(directIterator.Value()));
                }
            }

            std::vector<bool> matchedDirectUses(directUses.size(), false);
            int unexpectedExploredUses = 0;
            for (BRepTools_WireExplorer edgeExplorer(wire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                const TopoDS_Edge exploredEdge = edgeExplorer.Current();
                bool matched = false;
                for (std::size_t directIndex = 0;
                     directIndex < directUses.size();
                     ++directIndex) {
                    if (!matchedDirectUses[directIndex] &&
                        directUses[directIndex].IsSame(exploredEdge)) {
                        matchedDirectUses[directIndex] = true;
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    ++unexpectedExploredUses;
                }
            }

            traversalMismatch += unexpectedExploredUses;
            traversalMismatch += static_cast<int>(std::count(
                matchedDirectUses.begin(),
                matchedDirectUses.end(),
                false));
        }

        if (traversalMismatch > bestTraversalMismatch) {
            bestTraversalMismatch = traversalMismatch;
            plan.defectiveFace = relatedFace;
            sourceIsChainA = canReplaceChainA;
        }
    }

    if (plan.defectiveFace.IsNull() || bestTraversalMismatch <= 0) {
        errorMessage =
            "no candidate-related face has both an incomplete wire traversal and one complete source chain";
        return false;
    }
    plan.missingWireEdgeUses = bestTraversalMismatch;

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(currentShape, edgeFaceMap);

    auto findOtherFaces = [&](const TopoDS_Edge& edge) {
        std::vector<TopoDS_Face> otherFaces;
        if (edge.IsNull() || !edgeFaceMap.Contains(edge)) {
            return otherFaces;
        }
        const TopTools_ListOfShape& adjacentFaces = edgeFaceMap.FindFromKey(edge);
        for (TopTools_ListIteratorOfListOfShape faceIterator(adjacentFaces);
             faceIterator.More();
             faceIterator.Next()) {
            const TopoDS_Face adjacentFace = TopoDS::Face(faceIterator.Value());
            if (adjacentFace.IsNull() || adjacentFace.IsSame(plan.defectiveFace) ||
                TopologyUtils::IsSameFaceInVector(otherFaces, adjacentFace)) {
                continue;
            }
            otherFaces.push_back(adjacentFace);
        }
        return otherFaces;
    };

    std::vector<TopoDS_Edge> unorderedBoundaryEdges;
    auto addBoundaryEdge = [&](const TopoDS_Edge& edge) {
        if (!edge.IsNull() && !isSameEdgeInVector(unorderedBoundaryEdges, edge)) {
            unorderedBoundaryEdges.push_back(edge);
        }
    };
    auto addRemovedEdge = [&](const TopoDS_Edge& edge) {
        if (!edge.IsNull() && !isSameEdgeInVector(plan.removedEdges, edge)) {
            plan.removedEdges.push_back(edge);
        }
    };

    TopoDS_Face forwardDefectiveFace = plan.defectiveFace;
    forwardDefectiveFace.Orientation(TopAbs_FORWARD);
    std::vector<TopoDS_Edge> defectiveFaceEdges;
    for (TopExp_Explorer wireExplorer(forwardDefectiveFace, TopAbs_WIRE);
         wireExplorer.More();
         wireExplorer.Next()) {
        const TopoDS_Wire wire = TopoDS::Wire(wireExplorer.Current());
        for (TopoDS_Iterator edgeIterator(wire, Standard_False, Standard_False);
             edgeIterator.More();
             edgeIterator.Next()) {
            if (edgeIterator.Value().ShapeType() != TopAbs_EDGE) {
                continue;
            }
            const TopoDS_Edge edge = TopoDS::Edge(edgeIterator.Value());
            if (isSameEdgeInVector(defectiveFaceEdges, edge)) {
                errorMessage =
                    "defective face contains a repeated or seam edge use; delete-and-fill requires a simple boundary";
                return false;
            }
            defectiveFaceEdges.push_back(edge);

            if (!edgeFaceMap.Contains(edge)) {
                errorMessage = "defective face edge is missing from the current edge-face map";
                return false;
            }
            const std::vector<TopoDS_Face> otherFaces = findOtherFaces(edge);
            if (otherFaces.size() > 1) {
                errorMessage =
                    "deleting the defective face would retain a non-manifold boundary edge";
                return false;
            }
            if (otherFaces.empty()) {
                addRemovedEdge(edge);
            }
            else {
                addBoundaryEdge(edge);
            }
        }
    }

    const EdgeChain& sourceChain = sourceIsChainA
        ? candidate.chainA
        : candidate.chainB;
    const EdgeChain& oppositeChain = sourceIsChainA
        ? candidate.chainB
        : candidate.chainA;
    for (const TopoDS_Edge& sourceEdge : sourceChain.edges) {
        if (!isSameEdgeInVector(plan.removedEdges, sourceEdge)) {
            errorMessage =
                "source chain edge would survive defective face deletion";
            return false;
        }
    }

    // 缝隙另一侧的自由边在删除 face 后仍存在，并成为新补面的替代边界。
    for (const TopoDS_Edge& oppositeEdge : oppositeChain.edges) {
        const std::vector<TopoDS_Face> otherFaces = findOtherFaces(oppositeEdge);
        if (otherFaces.size() != 1) {
            errorMessage =
                "opposite chain edge does not have exactly one surviving owner face";
            return false;
        }
        addBoundaryEdge(oppositeEdge);
    }

    auto addSurvivingConnector = [&](const GapEnd& end) {
        if (end.type != GapEndType::ShortEdge || end.connectorEdge.IsNull()) {
            return true;
        }
        const std::vector<TopoDS_Face> otherFaces = findOtherFaces(end.connectorEdge);
        if (otherFaces.empty()) {
            return true;
        }
        if (otherFaces.size() != 1) {
            errorMessage =
                "short endpoint connector has ambiguous surviving owner faces";
            return false;
        }
        addBoundaryEdge(end.connectorEdge);
        return true;
    };
    if (!addSurvivingConnector(candidate.startEnd) ||
        !addSurvivingConnector(candidate.finishEnd)) {
        return false;
    }

    if (unorderedBoundaryEdges.size() < 3) {
        errorMessage = "delete-and-fill boundary has fewer than three edges";
        return false;
    }

    struct BoundaryLink {
        TopoDS_Edge edge;
        TopoDS_Vertex firstVertex;
        TopoDS_Vertex lastVertex;
        int firstIndex = 0;
        int lastIndex = 0;
    };
    TopTools_IndexedMapOfShape boundaryVertices;
    std::vector<BoundaryLink> links;
    links.reserve(unorderedBoundaryEdges.size());
    for (const TopoDS_Edge& edge : unorderedBoundaryEdges) {
        BoundaryLink link;
        link.edge = edge;
        if (!TopologyUtils::GetEdgeVertices(edge, link.firstVertex, link.lastVertex) ||
            link.firstVertex.IsSame(link.lastVertex)) {
            errorMessage =
                "delete-and-fill boundary contains an invalid or closed single edge";
            return false;
        }
        link.firstIndex = boundaryVertices.Add(link.firstVertex);
        link.lastIndex = boundaryVertices.Add(link.lastVertex);
        links.push_back(link);
    }

    std::vector<std::vector<std::size_t>> incidentLinks(
        static_cast<std::size_t>(boundaryVertices.Extent() + 1));
    for (std::size_t linkIndex = 0; linkIndex < links.size(); ++linkIndex) {
        incidentLinks[static_cast<std::size_t>(links[linkIndex].firstIndex)].push_back(linkIndex);
        incidentLinks[static_cast<std::size_t>(links[linkIndex].lastIndex)].push_back(linkIndex);
    }
    for (int vertexIndex = 1; vertexIndex <= boundaryVertices.Extent(); ++vertexIndex) {
        const std::size_t degree =
            incidentLinks[static_cast<std::size_t>(vertexIndex)].size();
        if (degree != 2) {
            std::ostringstream message;
            message << "surviving boundary is not one unbranched loop: vertexDegree="
                    << degree << ", expected=2";
            errorMessage = message.str();
            return false;
        }
    }

    std::vector<bool> usedLinks(links.size(), false);
    TopoDS_Vertex loopStart = links.front().firstVertex;
    TopoDS_Vertex currentVertex = loopStart;
    plan.boundary.vertices.push_back(loopStart);
    for (std::size_t step = 0; step < links.size(); ++step) {
        const int currentVertexIndex = boundaryVertices.FindIndex(currentVertex);
        if (currentVertexIndex <= 0) {
            errorMessage = "failed to locate the current boundary vertex";
            return false;
        }

        std::size_t nextLinkIndex = links.size();
        if (step == 0) {
            // 闭环起点天然有两个方向，固定从第一条边开始即可消除该二义性。
            nextLinkIndex = 0;
        }
        else {
            for (const std::size_t incidentLinkIndex :
                 incidentLinks[static_cast<std::size_t>(currentVertexIndex)]) {
                if (!usedLinks[incidentLinkIndex]) {
                    if (nextLinkIndex != links.size()) {
                        errorMessage =
                            "surviving boundary traversal has more than one next edge";
                        return false;
                    }
                    nextLinkIndex = incidentLinkIndex;
                }
            }
        }
        if (nextLinkIndex == links.size()) {
            errorMessage = "surviving boundary traversal ended before using every edge";
            return false;
        }

        const BoundaryLink& link = links[nextLinkIndex];
        TopoDS_Vertex nextVertex;
        if (link.firstVertex.IsSame(currentVertex)) {
            nextVertex = link.lastVertex;
        }
        else if (link.lastVertex.IsSame(currentVertex)) {
            nextVertex = link.firstVertex;
        }
        else {
            errorMessage = "selected boundary edge is not incident to the current vertex";
            return false;
        }

        plan.boundary.edges.push_back(
            TopologyUtils::OrientEdgeFromTo(link.edge, currentVertex, nextVertex));
        plan.boundary.vertices.push_back(nextVertex);
        plan.boundary.length += GeometryUtils::ComputeEdgeLength(link.edge);
        usedLinks[nextLinkIndex] = true;
        currentVertex = nextVertex;
    }

    if (!currentVertex.IsSame(loopStart) ||
        std::find(usedLinks.begin(), usedLinks.end(), false) != usedLinks.end()) {
        errorMessage =
            "surviving boundary contains multiple loops or is not topologically closed";
        return false;
    }

    for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
        const std::vector<TopoDS_Face> otherFaces = findOtherFaces(boundaryEdge);
        if (otherFaces.size() != 1) {
            errorMessage =
                "ordered boundary edge does not have exactly one surviving owner face";
            return false;
        }
        plan.boundary.ownerFaces.push_back(otherFaces.front());
    }
    return true;
}

bool HealClosedNarrowGap::buildDeleteAndFillPatchFace(
    const DeleteAndFillPlan& plan,
    TopoDS_Face& patchFace,
    double& g0Error,
    std::string& errorMessage) const
{
    patchFace.Nullify();
    g0Error = 0.0;
    if (plan.boundary.edges.size() < 3 ||
        plan.boundary.vertices.size() != plan.boundary.edges.size() + 1 ||
        !plan.boundary.vertices.front().IsSame(plan.boundary.vertices.back())) {
        errorMessage = "delete-and-fill boundary is not a valid closed chain";
        return false;
    }

    try {
        BRepFill_Filling filling;
        filling.SetResolParam(
            params_.deleteAndFillDegree,
            params_.deleteAndFillPointsOnCurve,
            params_.deleteAndFillIterations,
            Standard_False);
        filling.SetConstrParam(
            params_.deleteAndFillTolerance,
            params_.deleteAndFillTolerance,
            1.0e-2,
            1.0e-1);
        filling.SetApproxParam(
            params_.deleteAndFillMaxDegree,
            params_.deleteAndFillMaxSegments);

        // 直接以删除后留下的三维边界拟合新曲面，不使用旧 face 的 surface。
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            filling.Add(boundaryEdge, GeomAbs_C0, Standard_True);
        }
        filling.Build();
        if (!filling.IsDone()) {
            errorMessage = "BRepFill_Filling did not build a replacement face";
            return false;
        }

        const TopoDS_Face temporaryFillingFace = filling.Face();
        if (temporaryFillingFace.IsNull()) {
            errorMessage = "BRepFill_Filling returned a null replacement face";
            return false;
        }
        g0Error = filling.G0Error();
        if (!std::isfinite(g0Error) ||
            g0Error > params_.deleteAndFillMaxG0Error) {
            std::ostringstream message;
            message << "delete-and-fill boundary G0 error is too large: error="
                    << g0Error
                    << ", maximum=" << params_.deleteAndFillMaxG0Error;
            errorMessage = message.str();
            return false;
        }

        const Handle(Geom_Surface) locatedFillingSurface =
            BRep_Tool::Surface(temporaryFillingFace);
        if (locatedFillingSurface.IsNull()) {
            errorMessage = "delete-and-fill replacement face has no support surface";
            return false;
        }
        if (Handle(Geom_BSplineSurface)::DownCast(locatedFillingSurface).IsNull()) {
            errorMessage =
                "DeleteAndFillDefectiveFace did not produce a B-spline support surface";
            return false;
        }

        double pcurveTolerance = std::max(
            params_.deleteAndFillTolerance,
            std::max(g0Error, Precision::Confusion()));
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            pcurveTolerance = std::max(
                pcurveTolerance,
                BRep_Tool::Tolerance(boundaryEdge));
        }

        BRep_Builder builder;
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            Standard_Real firstParameter = 0.0;
            Standard_Real lastParameter = 0.0;
            Handle(Geom_Curve) curve3d =
                BRep_Tool::Curve(boundaryEdge, firstParameter, lastParameter);
            if (curve3d.IsNull() ||
                !std::isfinite(firstParameter) ||
                !std::isfinite(lastParameter) ||
                lastParameter - firstParameter <= Precision::PConfusion()) {
                errorMessage =
                    "surviving boundary edge has no valid 3d curve range";
                return false;
            }

            // 由原 3D 参数直接投影得到 same-parameter PCurve，避免 FixAddPCurve
            // 对新拟合 surface 只得到“非空但参数不一致”的二维曲线。
            ShapeConstruct_ProjectCurveOnSurface projector;
            projector.Init(locatedFillingSurface, pcurveTolerance);
            Handle(Geom2d_Curve) projectedPCurve;
            if (!projector.Perform(
                    curve3d,
                    firstParameter,
                    lastParameter,
                    projectedPCurve) &&
                !projector.PerformByProjLib(
                    curve3d,
                    firstParameter,
                    lastParameter,
                    projectedPCurve)) {
                errorMessage =
                    "failed to project a surviving boundary edge onto the filling surface";
                return false;
            }
            if (projectedPCurve.IsNull()) {
                errorMessage =
                    "boundary projection returned a null pcurve";
                return false;
            }

            builder.UpdateEdge(
                boundaryEdge,
                projectedPCurve,
                locatedFillingSurface,
                TopLoc_Location(),
                pcurveTolerance);
            builder.Range(
                boundaryEdge,
                locatedFillingSurface,
                TopLoc_Location(),
                firstParameter,
                lastParameter);
        }

        TopoDS_Wire boundaryWire;
        builder.MakeWire(boundaryWire);
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            builder.Add(boundaryWire, boundaryEdge);
        }
        boundaryWire.Closed(Standard_True);

        // Filling 只提供新 surface；最终 face 必须显式复用删除后留下的原 edge TShape。
        builder.MakeFace(
            patchFace,
            locatedFillingSurface,
            TopLoc_Location(),
            pcurveTolerance);
        patchFace.Orientation(TopAbs_FORWARD);
        builder.Add(patchFace, boundaryWire);
        builder.UpdateFace(patchFace, pcurveTolerance);

        // Filling 的 G0Error 针对其内部边界；最终复用原 edge 后必须重新测量真实偏差。
        double measuredTolerance = pcurveTolerance;
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            BRepLib_CheckCurveOnSurface checker(boundaryEdge, patchFace);
            checker.Perform();
            if (!checker.IsDone() || !std::isfinite(checker.MaxDistance())) {
                errorMessage =
                    "failed to measure delete-and-fill boundary deviation";
                return false;
            }

            const double existingEdgeTolerance =
                BRep_Tool::Tolerance(boundaryEdge);
            const double requiredEdgeTolerance = std::max(
                existingEdgeTolerance,
                checker.MaxDistance() + Precision::Confusion());
            const double allowedEdgeTolerance = std::max(
                existingEdgeTolerance,
                params_.deleteAndFillMaxG0Error);
            if (requiredEdgeTolerance > allowedEdgeTolerance) {
                std::ostringstream message;
                message << "rebuilt boundary deviation exceeds delete-and-fill tolerance: required="
                        << requiredEdgeTolerance
                        << ", maximum=" << allowedEdgeTolerance;
                errorMessage = message.str();
                return false;
            }

            builder.UpdateEdge(boundaryEdge, requiredEdgeTolerance);
            measuredTolerance = std::max(
                measuredTolerance,
                requiredEdgeTolerance);
        }
        pcurveTolerance = measuredTolerance;
        builder.UpdateFace(patchFace, pcurveTolerance);

        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            TopoDS_Edge boundaryUse = boundaryEdge;
            std::string pcurveError;
            if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                    boundaryUse,
                    patchFace,
                    pcurveTolerance,
                    &pcurveError)) {
                errorMessage =
                    "failed to prepare pcurve on delete-and-fill face: "
                    + pcurveError;
                return false;
            }
            if (!boundaryUse.IsSame(boundaryEdge)) {
                errorMessage =
                    "pcurve repair changed a surviving boundary edge identity";
                return false;
            }
        }

        TopTools_IndexedMapOfShape patchEdges;
        TopExp::MapShapes(patchFace, TopAbs_EDGE, patchEdges);
        if (patchEdges.Extent() != static_cast<int>(plan.boundary.edges.size())) {
            errorMessage =
                "rebuilt face boundary edge count differs from the surviving loop";
            return false;
        }
        for (const TopoDS_Edge& boundaryEdge : plan.boundary.edges) {
            bool edgeWasReused = false;
            for (int edgeIndex = 1; edgeIndex <= patchEdges.Extent(); ++edgeIndex) {
                if (patchEdges.FindKey(edgeIndex).IsSame(boundaryEdge)) {
                    edgeWasReused = true;
                    break;
                }
            }
            if (!edgeWasReused) {
                errorMessage =
                    "rebuilt face did not reuse every surviving boundary edge";
                return false;
            }
        }

        BRepCheck_Analyzer faceAnalyzer(patchFace, Standard_True);
        if (!faceAnalyzer.IsValid()) {
            // 同一闭环的另一侧可能相对拟合 surface 采用相反环向；只整体反转 wire，
            // 不改变任何边界 edge TShape，再重新检查 face。
            TopoDS_Wire reversedBoundaryWire = boundaryWire;
            reversedBoundaryWire.Reverse();
            TopoDS_Face reversedPatchFace;
            builder.MakeFace(
                reversedPatchFace,
                locatedFillingSurface,
                TopLoc_Location(),
                pcurveTolerance);
            reversedPatchFace.Orientation(TopAbs_FORWARD);
            builder.Add(reversedPatchFace, reversedBoundaryWire);
            builder.UpdateFace(reversedPatchFace, pcurveTolerance);

            BRepCheck_Analyzer reversedAnalyzer(
                reversedPatchFace,
                Standard_True);
            if (!reversedAnalyzer.IsValid()) {
                errorMessage = "delete-and-fill replacement face is invalid after trying both wire orientations: "
                    + RepairDiagnostics::DescribeBRepCheckFailures(
                        reversedAnalyzer,
                        reversedPatchFace);
                return false;
            }
            patchFace = reversedPatchFace;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage =
            "OCCT failure while building DeleteAndFillDefectiveFace patch: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::orientDeleteAndFillPatchFace(
    const DeleteAndFillPlan& plan,
    TopoDS_Face& patchFace,
    bool& reversedPatchFace,
    std::string& errorMessage) const
{
    reversedPatchFace = false;
    if (patchFace.IsNull() ||
        plan.boundary.edges.empty() ||
        plan.boundary.ownerFaces.size() != plan.boundary.edges.size()) {
        errorMessage =
            "delete-and-fill patch or surviving boundary data is invalid";
        return false;
    }

    // 先以支撑面的正向法线为基准比较。若每条共享边都与相邻 face 同向，
    // 说明整个 patch face 应翻转；混合结果则无法通过一次整体翻转修正。
    TopoDS_Face forwardPatchFace = patchFace;
    forwardPatchFace.Orientation(TopAbs_FORWARD);
    int sameDirectionCount = 0;
    int oppositeDirectionCount = 0;
    for (std::size_t edgeIndex = 0;
         edgeIndex < plan.boundary.edges.size();
         ++edgeIndex) {
        TopAbs_Orientation patchEdgeOrientation = TopAbs_FORWARD;
        TopAbs_Orientation ownerEdgeOrientation = TopAbs_FORWARD;
        std::string orientationError;
        if (!findUniqueEdgeUseOrientation(
                forwardPatchFace,
                plan.boundary.edges[edgeIndex],
                patchEdgeOrientation,
                orientationError)) {
            errorMessage = "failed to read patch boundary orientation: "
                + orientationError;
            return false;
        }
        if (!findUniqueEdgeUseOrientation(
                plan.boundary.ownerFaces[edgeIndex],
                plan.boundary.edges[edgeIndex],
                ownerEdgeOrientation,
                orientationError)) {
            errorMessage = "failed to read surviving face boundary orientation: "
                + orientationError;
            return false;
        }

        if (patchEdgeOrientation == ownerEdgeOrientation) {
            ++sameDirectionCount;
        }
        else {
            ++oppositeDirectionCount;
        }
    }

    if (sameDirectionCount > 0 && oppositeDirectionCount > 0) {
        std::ostringstream message;
        message << "surviving faces require inconsistent patch orientations: sameDirectionEdges="
                << sameDirectionCount
                << ", oppositeDirectionEdges=" << oppositeDirectionCount;
        errorMessage = message.str();
        return false;
    }

    reversedPatchFace = sameDirectionCount > 0;
    patchFace.Orientation(
        reversedPatchFace ? TopAbs_REVERSED : TopAbs_FORWARD);
    return true;
}

bool HealClosedNarrowGap::validateDeleteAndFillBoundaryOrientations(
    const TopoDS_Shape& shape,
    const EdgeChain& boundary,
    std::string& errorMessage) const
{
    if (shape.IsNull() || boundary.edges.empty()) {
        errorMessage =
            "shape or delete-and-fill boundary is invalid for orientation validation";
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);
    for (std::size_t edgeIndex = 0;
         edgeIndex < boundary.edges.size();
         ++edgeIndex) {
        const TopoDS_Edge& boundaryEdge = boundary.edges[edgeIndex];
        if (boundaryEdge.IsNull() || !edgeFaceMap.Contains(boundaryEdge)) {
            errorMessage =
                "delete-and-fill boundary edge is missing during orientation validation";
            return false;
        }

        const TopTools_ListOfShape& adjacentFaces =
            edgeFaceMap.FindFromKey(boundaryEdge);
        if (adjacentFaces.Extent() != 2) {
            std::ostringstream message;
            message << "delete-and-fill boundary edge does not have exactly two adjacent faces: edgeIndex="
                    << edgeIndex
                    << ", adjacentFaces=" << adjacentFaces.Extent();
            errorMessage = message.str();
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
            std::string orientationError;
            if (!findUniqueEdgeUseOrientation(
                    TopoDS::Face(faceIterator.Value()),
                    boundaryEdge,
                    orientations[faceIndex],
                    orientationError)) {
                errorMessage =
                    "failed to validate delete-and-fill boundary orientation: "
                    + orientationError;
                return false;
            }
        }

        if (orientations[0] == orientations[1]) {
            std::ostringstream message;
            message << "delete-and-fill boundary edge has the same orientation in both adjacent faces: edgeIndex="
                    << edgeIndex;
            errorMessage = message.str();
            return false;
        }
    }
    return true;
}

bool HealClosedNarrowGap::tryRepairByDeletingAndFillingDefectiveFace(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const ClosedGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    DeleteAndFillPlan plan;
    if (!buildDeleteAndFillPlan(
            currentShape,
            candidate,
            plan,
            errorMessage)) {
        return false;
    }

    {
        std::ostringstream message;
        message << "DeleteAndFillDefectiveFace plan: missingWireEdgeUses="
                << plan.missingWireEdgeUses
                << ", boundaryEdges=" << plan.boundary.edges.size()
                << ", removedEdges=" << plan.removedEdges.size();
        appendMessage(message.str());
    }

    TopoDS_Face patchFace;
    double g0Error = 0.0;
    if (!buildDeleteAndFillPatchFace(
            plan,
            patchFace,
            g0Error,
            errorMessage)) {
        return false;
    }

    bool reversedPatchFace = false;
    if (!orientDeleteAndFillPatchFace(
            plan,
            patchFace,
            reversedPatchFace,
            errorMessage)) {
        errorMessage = "failed to orient DeleteAndFill patch face: "
            + errorMessage;
        return false;
    }
    {
        std::ostringstream message;
        message << "DeleteAndFill patch orientation: reversedAgainstSurface="
                << (reversedPatchFace ? "true" : "false")
                << ", boundaryEdges=" << plan.boundary.edges.size();
        appendMessage(message.str());
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    TopoDS_Shape trialShape;
    TopoDS_Shape locallyReplacedShape;
    try {
        // source key 仍规范为 FORWARD，但 replacement value 编码“旧 face use ->
        // 已按 surviving faces 定向的 patch face use”的相对方向，不能继续盲目继承旧方向。
        TopoDS_Face sourceKey = plan.defectiveFace;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Face replacementValue = patchFace;
        const TopAbs_Orientation sourceUseOrientation =
            plan.defectiveFace.Orientation();
        const TopAbs_Orientation desiredPatchOrientation =
            patchFace.Orientation();
        if ((sourceUseOrientation != TopAbs_FORWARD &&
             sourceUseOrientation != TopAbs_REVERSED) ||
            (desiredPatchOrientation != TopAbs_FORWARD &&
             desiredPatchOrientation != TopAbs_REVERSED)) {
            errorMessage =
                "defective or replacement face has unsupported orientation";
            return false;
        }
        replacementValue.Orientation(
            sourceUseOrientation == desiredPatchOrientation
                ? TopAbs_FORWARD
                : TopAbs_REVERSED);

        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage =
                "replacing the defective face produced a null shape";
            return false;
        }

        std::string boundaryOrientationError;
        if (!validateDeleteAndFillBoundaryOrientations(
                trialShape,
                plan.boundary,
                boundaryOrientationError)) {
            errorMessage =
                "DeleteAndFill local replacement has inconsistent face orientation: "
                + boundaryOrientationError;
            return false;
        }
        locallyReplacedShape = trialShape;

        appendMessage(
            "DeleteAndFill adjacency after local face replacement: "
            + RepairDiagnostics::DescribeEdgeFaceCounts(
                locallyReplacedShape,
                plan.boundary.edges));

        if (params_.runShapeFixAfterRepair) {
            ShapeFix_Shape fixer(trialShape);
            const double fixPrecision = std::max(
                Precision::Confusion(),
                params_.deleteAndFillTolerance);
            fixer.SetPrecision(fixPrecision);
            fixer.SetMaxTolerance(std::max(
                fixPrecision,
                params_.deleteAndFillMaxG0Error));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            const TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage =
            "OCCT failure while replacing the defective face: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    appendMessage(
        std::string(params_.runShapeFixAfterRepair
            ? "DeleteAndFill adjacency after final ShapeFix: "
            : "DeleteAndFill adjacency before validation: ")
        + RepairDiagnostics::DescribeEdgeFaceCounts(trialShape, plan.boundary.edges));

    auto defectiveFaceWasRemoved = [&](const TopoDS_Shape& shape) {
        TopTools_IndexedMapOfShape resultFaces;
        TopExp::MapShapes(shape, TopAbs_FACE, resultFaces);
        for (int faceIndex = 1; faceIndex <= resultFaces.Extent(); ++faceIndex) {
            if (resultFaces.FindKey(faceIndex).IsSame(plan.defectiveFace)) {
                return false;
            }
        }
        return true;
    };
    if (!defectiveFaceWasRemoved(trialShape)) {
        errorMessage = "defective face is still present after replacement";
        return false;
    }

    auto validateDeleteAndFillResult = [&](
        const TopoDS_Shape& shape,
        std::string& validationError) {
        if (!validateRepairResult(
                existingFailures,
                shape,
                plan.boundary,
                plan.removedEdges,
                freeEdgesBefore,
                validationError)) {
            return false;
        }

        std::string boundaryOrientationError;
        if (!validateDeleteAndFillBoundaryOrientations(
                shape,
                plan.boundary,
                boundaryOrientationError)) {
            validationError =
                "DeleteAndFill boundary orientation validation failed: "
                + boundaryOrientationError;
            return false;
        }
        return true;
    };

    std::string validationError;
    if (!validateDeleteAndFillResult(trialShape, validationError)) {
        std::string localValidationError;
        const bool canUseLocalResult =
            params_.runShapeFixAfterRepair &&
            !locallyReplacedShape.IsNull() &&
            !trialShape.IsSame(locallyReplacedShape) &&
            defectiveFaceWasRemoved(locallyReplacedShape) &&
            validateDeleteAndFillResult(
                locallyReplacedShape,
                localValidationError);
        if (!canUseLocalResult) {
            errorMessage = validationError;
            if (!localValidationError.empty()) {
                errorMessage += "\nlocal replacement validation also failed: "
                    + localValidationError;
            }
            return false;
        }

        appendMessage(
            "Warning: final ShapeFix result was rejected; using the validated local DeleteAndFill result: "
            + validationError);
        trialShape = locallyReplacedShape;
    }

    {
        std::ostringstream message;
        message << "Candidate repair accepted by DeleteAndFillDefectiveFace: boundaryEdges="
                << plan.boundary.edges.size()
                << ", removedEdges=" << plan.removedEdges.size()
                << ", g0Error=" << g0Error
                << ", freeEdgesBefore=" << freeEdgesBefore
                << ", freeEdgesAfter="
                << TopologyUtils::CountFreeEdges(trialShape);
        appendMessage(message.str());
    }
    repairedShape = trialShape;
    return true;
}

