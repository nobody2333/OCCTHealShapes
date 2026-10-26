#include "HealClosedNarrowGap.h"
#include "RepairDiagnostics.h"

#include "TopologyUtils.h"

#include <algorithm>
#include <limits>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>

bool HealClosedNarrowGap::rebuildFacesByCollapsingToMidline(
    const TopoDS_Shape& currentShape,
    const ClosedGapCandidate& candidate,
    const MidlineChainData& midlineData,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }
    if (midlineData.chain.edges.empty()) {
        errorMessage = "midline chain is empty";
        return false;
    }
    if (midlineData.chain.vertices.size() !=
        midlineData.chain.edges.size() + 1) {
        errorMessage = "midline target vertex and edge counts are inconsistent";
        return false;
    }
    if (!validateReplacementMapping(
            midlineData.replacementsForChainA,
            candidate.chainA.edges.size(),
            midlineData.chain.edges.size(),
            "chainA",
            errorMessage) ||
        !validateReplacementMapping(
            midlineData.replacementsForChainB,
            candidate.chainB.edges.size(),
            midlineData.chain.edges.size(),
            "chainB",
            errorMessage)) {
        return false;
    }

    std::vector<FaceReplacementPlan> facePlans;
    auto findOrAddFacePlan = [&](const TopoDS_Face& face) -> FaceReplacementPlan* {
        for (FaceReplacementPlan& plan : facePlans) {
            if (plan.face.IsSame(face)) {
                return &plan;
            }
        }

        FaceReplacementPlan plan;
        plan.face = face;
        facePlans.push_back(plan);
        return &facePlans.back();
    };

    // 按原 edge 的 owner face 登记局部替换，不再要求完整 chain 位于单个 face/wire 中。
    auto addChainTasks = [&](const EdgeChain& sourceChain,
                             const std::vector<std::vector<std::size_t>>& replacements,
                             const char* chainName) {
        if (sourceChain.ownerFaces.size() != sourceChain.edges.size() ||
            sourceChain.vertices.size() != sourceChain.edges.size() + 1 ||
            replacements.size() != sourceChain.edges.size()) {
            errorMessage = std::string(chainName) + " replacement mapping is inconsistent";
            return false;
        }

        for (std::size_t edgeIndex = 0; edgeIndex < sourceChain.edges.size(); ++edgeIndex) {
            const TopoDS_Face& ownerFace = sourceChain.ownerFaces[edgeIndex];
            if (ownerFace.IsNull()) {
                errorMessage = std::string(chainName) + " contains an edge without owner face";
                return false;
            }
            if (replacements[edgeIndex].empty()) {
                errorMessage = std::string(chainName) + " edge has no midline replacement segment";
                return false;
            }
            for (std::size_t replacementIndex = 0;
                 replacementIndex < replacements[edgeIndex].size();
                 ++replacementIndex) {
                const std::size_t targetEdgeIndex =
                    replacements[edgeIndex][replacementIndex];
                if (targetEdgeIndex >= midlineData.chain.edges.size()) {
                    errorMessage = std::string(chainName)
                        + " replacement target edge index is out of range";
                    return false;
                }
                if (replacementIndex > 0 &&
                    targetEdgeIndex !=
                        replacements[edgeIndex][replacementIndex - 1] + 1) {
                    errorMessage = std::string(chainName)
                        + " replacement target edge indices are not continuous";
                    return false;
                }
            }

            FaceReplacementPlan* plan = findOrAddFacePlan(ownerFace);
            for (const EdgeReplacementTask& existingTask : plan->edgeTasks) {
                if (existingTask.sourceEdge.IsSame(sourceChain.edges[edgeIndex])) {
                    errorMessage = "the same source edge appears more than once in midline replacement tasks";
                    return false;
                }
                for (const std::size_t existingTargetIndex :
                     existingTask.replacementEdgeIndices) {
                    if (std::find(
                            replacements[edgeIndex].begin(),
                            replacements[edgeIndex].end(),
                            existingTargetIndex)
                        != replacements[edgeIndex].end()) {
                        errorMessage = "the same canonical target edge would be inserted "
                            "more than once into one face";
                        return false;
                    }
                }
            }

            EdgeReplacementTask task;
            task.sourceEdge = sourceChain.edges[edgeIndex];
            task.sourceStart = sourceChain.vertices[edgeIndex];
            task.sourceEnd = sourceChain.vertices[edgeIndex + 1];
            task.replacementEdgeIndices = replacements[edgeIndex];
            plan->edgeTasks.push_back(task);
        }
        return true;
    };

    if (!addChainTasks(
            candidate.chainA,
            midlineData.replacementsForChainA,
            "chainA") ||
        !addChainTasks(
            candidate.chainB,
            midlineData.replacementsForChainB,
            "chainB")) {
        return false;
    }

    auto addEndEdgeRemoval = [&](const GapEnd& end) {
        if (end.type == GapEndType::SharedTip) {
            return true;
        }
        if (end.type != GapEndType::ShortEdge ||
            end.connectorEdge.IsNull() ||
            end.connectorOwnerFace.IsNull()) {
            errorMessage = "short endpoint removal data is incomplete";
            return false;
        }

        FaceReplacementPlan* plan =
            findOrAddFacePlan(end.connectorOwnerFace);
        for (const EdgeReplacementTask& edgeTask : plan->edgeTasks) {
            if (edgeTask.sourceEdge.IsSame(end.connectorEdge)) {
                errorMessage = "short endpoint edge is also registered as a midline replacement edge";
                return false;
            }
        }
        for (const TopoDS_Edge& existingEdge : plan->edgesToRemove) {
            if (existingEdge.IsSame(end.connectorEdge)) {
                errorMessage = "the same short endpoint edge is registered more than once";
                return false;
            }
        }
        plan->edgesToRemove.push_back(end.connectorEdge);
        return true;
    };

    if (!addEndEdgeRemoval(candidate.startEnd) ||
        !addEndEdgeRemoval(candidate.finishEnd)) {
        return false;
    }

    if (facePlans.empty()) {
        errorMessage = "closed-gap chains have no owner faces";
        return false;
    }

    try {
        BRepTools_ReShape topologyReshaper;
        BRepTools_ReShape vertexReshaper;
        std::vector<TopoDS_Vertex> replacedSourceVertices;
        TopTools_MapOfShape uniqueSourceVertices;
        TopTools_MapOfShape endpointEdgesToRemove;
        for (const FaceReplacementPlan& plan : facePlans) {
            for (const TopoDS_Edge& edgeToRemove : plan.edgesToRemove) {
                endpointEdgesToRemove.Add(edgeToRemove);
            }
        }

        auto registerVertexReplacement = [&](const TopoDS_Vertex& sourceVertex,
                                             const TopoDS_Vertex& midlineVertex) {
            if (sourceVertex.IsNull() || midlineVertex.IsNull()) {
                errorMessage = "cannot register a null midline vertex replacement";
                return false;
            }
            if (sourceVertex.IsSame(midlineVertex)) {
                return true;
            }

            // ReShape 的替换键和值必须去除 edge use 携带的局部方向，
            // 否则 REVERSED source vertex 会使 replacement use 的方向整体翻转。
            TopoDS_Vertex normalizedSourceVertex = sourceVertex;
            normalizedSourceVertex.Orientation(TopAbs_FORWARD);
            TopoDS_Vertex normalizedMidlineVertex = midlineVertex;
            normalizedMidlineVertex.Orientation(TopAbs_FORWARD);

            if (vertexReshaper.IsRecorded(normalizedSourceVertex)) {
                const TopoDS_Shape existingReplacement =
                    vertexReshaper.Value(normalizedSourceVertex);
                if (existingReplacement.IsNull() ||
                    !existingReplacement.IsSame(normalizedMidlineVertex)) {
                    errorMessage = "one source vertex maps to different midline vertices";
                    return false;
                }
                return true;
            }

            // 两侧对应旧顶点统一替换成同一个规范中点 vertex。
            vertexReshaper.Replace(normalizedSourceVertex, normalizedMidlineVertex);
            topologyReshaper.Replace(normalizedSourceVertex, normalizedMidlineVertex);
            if (uniqueSourceVertices.Add(normalizedSourceVertex)) {
                replacedSourceVertices.push_back(normalizedSourceVertex);
            }
            return true;
        };

        const std::size_t invalidTargetVertexIndex =
            std::numeric_limits<std::size_t>::max();
        auto registerChainVertexReplacements =
            [&](const EdgeChain& sourceChain,
                const std::vector<std::size_t>& targetVertexIndices,
                const char* chainName) {
                if (targetVertexIndices.size() !=
                    sourceChain.vertices.size()) {
                    errorMessage = std::string(
                        "vertex replacement failure: ")
                        + chainName
                        + " source vertex mapping count is inconsistent";
                    return false;
                }
                if (targetVertexIndices.empty() ||
                    targetVertexIndices.front() != 0 ||
                    targetVertexIndices.back() ==
                        invalidTargetVertexIndex ||
                    targetVertexIndices.back() !=
                        midlineData.chain.vertices.size() - 1) {
                    errorMessage = std::string(
                        "vertex replacement failure: ")
                        + chainName
                        + " endpoint mapping is inconsistent";
                    return false;
                }

                // A、B 两侧分别按各自映射替换。共享尖点会因 IsSame 直接跳过，
                // 短端边的两个端点则会映射到同一个中点 vertex。
                for (std::size_t sourceVertexIndex = 0;
                     sourceVertexIndex < sourceChain.vertices.size();
                     ++sourceVertexIndex) {
                    const std::size_t targetVertexIndex =
                        targetVertexIndices[sourceVertexIndex];
                    if (targetVertexIndex ==
                            invalidTargetVertexIndex ||
                        targetVertexIndex >=
                            midlineData.chain.vertices.size() ||
                        (sourceVertexIndex == 0 &&
                            targetVertexIndex != 0) ||
                        (sourceVertexIndex + 1 ==
                                sourceChain.vertices.size() &&
                            targetVertexIndex + 1 !=
                                midlineData.chain.vertices.size()) ||
                        (sourceVertexIndex > 0 &&
                            sourceVertexIndex + 1 <
                                sourceChain.vertices.size() &&
                            (targetVertexIndex == 0 ||
                                targetVertexIndex + 1 >=
                                    midlineData.chain.vertices.size()))) {
                        errorMessage = std::string(
                            "vertex replacement failure: ")
                            + chainName
                            + " target vertex index is invalid";
                        return false;
                    }

                    if (!registerVertexReplacement(
                            sourceChain.vertices[sourceVertexIndex],
                            midlineData.chain.vertices[
                                targetVertexIndex])) {
                        if (errorMessage.rfind(
                                "vertex replacement failure:",
                                0) != 0) {
                            errorMessage =
                                "vertex replacement failure: "
                                + errorMessage;
                        }
                        return false;
                    }
                }
                return true;
            };

        if (!registerChainVertexReplacements(
                candidate.chainA,
                midlineData.targetVertexIndicesForChainA,
                "chainA") ||
            !registerChainVertexReplacements(
                candidate.chainB,
                midlineData.targetVertexIndicesForChainB,
                "chainB")) {
            return false;
        }

        // 先为所有连接待替换旧顶点的 edge 创建唯一替代品，确保相邻 face 继续共享同一 edge TShape。
        TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap;
        TopologyUtils::BuildVertexEdgeMap(currentShape, vertexEdgeMap);
        TopTools_MapOfShape transformedEdges;
        for (const TopoDS_Vertex& sourceVertex : replacedSourceVertices) {
            if (!vertexEdgeMap.Contains(sourceVertex)) {
                continue;
            }

            const TopTools_ListOfShape& incidentEdges = vertexEdgeMap.FindFromKey(sourceVertex);
            for (TopTools_ListIteratorOfListOfShape edgeIt(incidentEdges);
                 edgeIt.More();
                 edgeIt.Next()) {
                const TopoDS_Edge sourceEdge = TopoDS::Edge(edgeIt.Value());
                if (sourceEdge.IsNull() || !transformedEdges.Add(sourceEdge)) {
                    continue;
                }
                if (endpointEdgesToRemove.Contains(sourceEdge)) {
                    continue;
                }

                const TopoDS_Shape transformedShape = vertexReshaper.Apply(sourceEdge);
                if (transformedShape.IsNull() || transformedShape.ShapeType() != TopAbs_EDGE) {
                    errorMessage = "failed to propagate a midpoint vertex to an incident edge";
                    return false;
                }
                TopoDS_Edge transformedEdge = TopoDS::Edge(transformedShape);
                transformedEdge.Orientation(sourceEdge.Orientation());
                topologyReshaper.Replace(sourceEdge, transformedEdge);
            }
        }

        // 短端边收缩为一个 vertex 后不能保留为零长度 edge，必须显式删除。
        TopTools_MapOfShape removedEndpointEdges;
        for (const FaceReplacementPlan& plan : facePlans) {
            for (const TopoDS_Edge& edgeToRemove : plan.edgesToRemove) {
                if (removedEndpointEdges.Add(edgeToRemove)) {
                    topologyReshaper.Remove(edgeToRemove);
                }
            }
        }

        const double tolerance = std::max(params_.midlineRepairTolerance, Precision::Confusion());
        for (const FaceReplacementPlan& plan : facePlans) {
            TopoDS_Face rebuiltFace;
            if (!rebuildFaceByReplacingMappedEdges(
                    plan.face,
                    plan.edgeTasks,
                    plan.edgesToRemove,
                    midlineData.chain,
                    topologyReshaper,
                    tolerance,
                    rebuiltFace,
                    errorMessage)) {
                return false;
            }

            topologyReshaper.Replace(plan.face, rebuiltFace);
        }

        repairedShape = topologyReshaper.Apply(currentShape);
        if (repairedShape.IsNull()) {
            errorMessage = "BRepTools_ReShape produced a null shape after midline collapse";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "midline shape replacement failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFaceByReplacingMappedEdges(
    const TopoDS_Face& faceToRebuild,
    const std::vector<EdgeReplacementTask>& edgeTasks,
    const std::vector<TopoDS_Edge>& edgesToRemove,
    const EdgeChain& midlineChain,
    BRepTools_ReShape& topologyReshaper,
    double tolerance,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage)
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull()) {
        errorMessage = "face to rebuild is null";
        return false;
    }
    if (edgeTasks.empty() && edgesToRemove.empty()) {
        errorMessage = "face has no midline replacement or endpoint removal tasks";
        return false;
    }
    if (midlineChain.edges.empty()) {
        errorMessage = "canonical midline chain is empty";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "face to rebuild has no support surface";
            return false;
        }

        const double faceTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(faceToRebuild),
            tolerance
        });

        BRep_Builder builder;
        TopoDS_Face rebuilt;
        builder.MakeFace(rebuilt, surface, surfaceLocation, faceTolerance);
        builder.NaturalRestriction(rebuilt, BRep_Tool::NaturalRestriction(forwardFace));
        rebuilt.Orientation(TopAbs_FORWARD);

        std::vector<bool> replacedTasks(edgeTasks.size(), false);
        std::vector<bool> removedTasks(edgesToRemove.size(), false);
        auto findTaskIndex = [&](const TopoDS_Edge& edge) {
            for (std::size_t taskIndex = 0; taskIndex < edgeTasks.size(); ++taskIndex) {
                if (edgeTasks[taskIndex].sourceEdge.IsSame(edge)) {
                    return static_cast<int>(taskIndex);
                }
            }
            return -1;
        };
        auto findRemovalIndex = [&](const TopoDS_Edge& edge) {
            for (std::size_t removalIndex = 0;
                 removalIndex < edgesToRemove.size();
                 ++removalIndex) {
                if (edgesToRemove[removalIndex].IsSame(edge)) {
                    return static_cast<int>(removalIndex);
                }
            }
            return -1;
        };
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            // 直接子 edge use 是旧 wire 的完整基准，不能只依赖 WireExplorer 的返回结果。
            std::vector<TopoDS_Edge> directEdgeUses;
            bool directWireHasModification = false;
            for (TopoDS_Iterator directIterator(oldWire, Standard_False, Standard_False);
                 directIterator.More();
                 directIterator.Next()) {
                if (directIterator.Value().ShapeType() != TopAbs_EDGE) {
                    continue;
                }

                const TopoDS_Edge directEdge = TopoDS::Edge(directIterator.Value());
                directEdgeUses.push_back(directEdge);
                directWireHasModification =
                    directWireHasModification ||
                    findTaskIndex(directEdge) >= 0 ||
                    findRemovalIndex(directEdge) >= 0;
            }

            std::vector<TopoDS_Edge> wireEdges;
            bool wireHasModification = false;
            for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                const TopoDS_Edge edgeInWire = edgeExplorer.Current();
                wireEdges.push_back(edgeInWire);
                wireHasModification =
                    wireHasModification ||
                    findTaskIndex(edgeInWire) >= 0 ||
                    findRemovalIndex(edgeInWire) >= 0;
            }

            if (directWireHasModification) {
                // 按 IsSame 做多重集匹配，既检查数量，也检查 WireExplorer 是否遗漏具体 edge use。
                std::vector<bool> matchedDirectUses(directEdgeUses.size(), false);
                int unexpectedExploredUses = 0;
                for (const TopoDS_Edge& exploredEdge : wireEdges) {
                    bool matched = false;
                    for (std::size_t directIndex = 0;
                         directIndex < directEdgeUses.size();
                         ++directIndex) {
                        if (!matchedDirectUses[directIndex] &&
                            directEdgeUses[directIndex].IsSame(exploredEdge)) {
                            matchedDirectUses[directIndex] = true;
                            matched = true;
                            break;
                        }
                    }
                    if (!matched) {
                        ++unexpectedExploredUses;
                    }
                }

                const int missingDirectUses = static_cast<int>(std::count(
                    matchedDirectUses.begin(), matchedDirectUses.end(), false));
                {
                    std::ostringstream os;
                    os << "Midline wire traversal check: directEdgeUses="
                       << directEdgeUses.size()
                       << ", exploredEdgeUses=" << wireEdges.size()
                       << ", missingDirectUses=" << missingDirectUses
                       << ", unexpectedExploredUses=" << unexpectedExploredUses;
                    appendMessage(os.str());
                }

                if (missingDirectUses > 0 || unexpectedExploredUses > 0 ||
                    !wireHasModification) {
                    std::ostringstream os;
                    os << "BRepTools_WireExplorer did not cover every direct edge use of the affected wire: "
                       << "direct=" << directEdgeUses.size()
                       << ", explored=" << wireEdges.size()
                       << ", missing=" << missingDirectUses
                       << ", unexpected=" << unexpectedExploredUses;
                    errorMessage = os.str();
                    return false;
                }
            }

            TopoDS_Wire wireToAdd = oldWire;
            if (wireHasModification) {
                std::vector<TopoDS_Edge> rebuiltEdgeUses;
                std::vector<TopoDS_Edge> expectedReplacementUses;
                for (const TopoDS_Edge& edgeInWire : wireEdges) {
                    const int removalIndex = findRemovalIndex(edgeInWire);
                    if (removalIndex >= 0) {
                        const std::size_t mappedRemovalIndex =
                            static_cast<std::size_t>(removalIndex);
                        if (removedTasks[mappedRemovalIndex]) {
                            errorMessage = "a short endpoint edge appears more than once in face wires";
                            return false;
                        }

                        // 两端点已经映射到同一个中点；短连接 edge 在新 wire 中直接省略。
                        removedTasks[mappedRemovalIndex] = true;
                        continue;
                    }

                    const int taskIndex = findTaskIndex(edgeInWire);
                    if (taskIndex < 0) {
                        const TopoDS_Shape adjustedShape = topologyReshaper.Apply(edgeInWire);
                        if (adjustedShape.IsNull() || adjustedShape.ShapeType() != TopAbs_EDGE) {
                            errorMessage = "failed to apply midpoint vertices to a neighboring wire edge";
                            return false;
                        }
                        TopoDS_Edge adjustedEdge = TopoDS::Edge(adjustedShape);
                        adjustedEdge.Orientation(edgeInWire.Orientation());
                        rebuiltEdgeUses.push_back(adjustedEdge);
                        continue;
                    }

                    const std::size_t mappedIndex = static_cast<std::size_t>(taskIndex);
                    if (replacedTasks[mappedIndex]) {
                        errorMessage = "a source edge appears more than once in face wires";
                        return false;
                    }

                    const EdgeReplacementTask& task = edgeTasks[mappedIndex];
                    TopoDS_Vertex wireStart;
                    TopoDS_Vertex wireEnd;
                    if (!TopologyUtils::GetEdgeVertices(edgeInWire, wireStart, wireEnd)) {
                        errorMessage = "failed to obtain source edge vertices from face wire";
                        return false;
                    }

                    bool reverseReplacement = false;
                    if (wireStart.IsSame(task.sourceStart) && wireEnd.IsSame(task.sourceEnd)) {
                        reverseReplacement = false;
                    }
                    else if (wireStart.IsSame(task.sourceEnd) && wireEnd.IsSame(task.sourceStart)) {
                        reverseReplacement = true;
                    }
                    else {
                        errorMessage = "source edge vertices do not match its chain interval";
                        return false;
                    }

                    auto addReplacementUse = [&](std::size_t replacementEdgeIndex,
                                                 bool reverseUse) {
                        if (replacementEdgeIndex >= midlineChain.edges.size()) {
                            errorMessage = "mapped midline edge index is out of range";
                            return false;
                        }

                        // 所有 face 都从规范化中线 chain 取同一个底层 edge，只调整当前 wire use 的方向。
                        TopoDS_Edge replacementUse =
                            midlineChain.edges[replacementEdgeIndex];
                        if (reverseUse) {
                            replacementUse.Reverse();
                        }

                        // 先为当前 face 补建并验证 pcurve，再把同一个底层中线 edge use 加入 wire。
                        if (!ensurePCurveOnFace(
                                replacementUse,
                                forwardFace,
                                faceTolerance,
                                errorMessage)) {
                            errorMessage = "failed to repair mapped midline pcurve: " + errorMessage;
                            return false;
                        }
                        if (!replacementUse.IsSame(
                                midlineChain.edges[replacementEdgeIndex])) {
                            errorMessage =
                                "pcurve repair changed the canonical midline edge identity";
                            return false;
                        }

                        rebuiltEdgeUses.push_back(replacementUse);
                        expectedReplacementUses.push_back(replacementUse);
                        return true;
                    };

                    if (reverseReplacement) {
                        // source edge 在 wire 中反向使用时，规范目标链的遍历方向也整体反转，
                        // 因此既要倒序插入 replacement edges，也要反转每个 edge use。
                        for (int replacementIndex =
                                 static_cast<int>(task.replacementEdgeIndices.size()) - 1;
                             replacementIndex >= 0;
                             --replacementIndex) {
                            if (!addReplacementUse(
                                    task.replacementEdgeIndices[
                                        static_cast<std::size_t>(replacementIndex)],
                                    true)) {
                                return false;
                            }
                        }
                    }
                    else {
                        for (const std::size_t replacementEdgeIndex :
                             task.replacementEdgeIndices) {
                            if (!addReplacementUse(replacementEdgeIndex, false)) {
                                return false;
                            }
                        }
                    }

                    replacedTasks[mappedIndex] = true;
                }

                if (rebuiltEdgeUses.empty()) {
                    errorMessage =
                        "affected wire is empty after mapped midline replacement";
                    return false;
                }

                // MakeWire 会在顶点仅几何接近时复制 edge 并替换端点，这会破坏两侧 face
                // 对 canonical midline edge TShape 的共享。这里先要求真实拓扑连接，再直接组装 wire。
                auto describeVertexConnection = [](const TopoDS_Vertex& lhs,
                                                   const TopoDS_Vertex& rhs) {
                    std::ostringstream os;
                    try {
                        os << "distance=" << BRep_Tool::Pnt(lhs).Distance(BRep_Tool::Pnt(rhs))
                           << ", lhsTolerance=" << BRep_Tool::Tolerance(lhs)
                           << ", rhsTolerance=" << BRep_Tool::Tolerance(rhs);
                    }
                    catch (const Standard_Failure&) {
                        os << "distance unavailable";
                    }
                    return os.str();
                };

                TopoDS_Vertex firstWireVertex;
                TopoDS_Vertex previousEndVertex;
                for (std::size_t edgeUseIndex = 0;
                     edgeUseIndex < rebuiltEdgeUses.size();
                     ++edgeUseIndex) {
                    TopoDS_Vertex edgeStart;
                    TopoDS_Vertex edgeEnd;
                    if (!TopologyUtils::GetEdgeVertices(
                            rebuiltEdgeUses[edgeUseIndex],
                            edgeStart,
                            edgeEnd)) {
                        errorMessage =
                            "failed to obtain rebuilt wire edge vertices";
                        return false;
                    }

                    if (edgeUseIndex == 0) {
                        firstWireVertex = edgeStart;
                    }
                    else if (!previousEndVertex.IsSame(edgeStart)) {
                        errorMessage =
                            "rebuilt wire edges are only geometrically connected, not topologically shared: "
                            + describeVertexConnection(previousEndVertex, edgeStart);
                        return false;
                    }
                    previousEndVertex = edgeEnd;
                }

                if (firstWireVertex.IsNull() || previousEndVertex.IsNull() ||
                    !previousEndVertex.IsSame(firstWireVertex)) {
                    errorMessage =
                        "rebuilt wire is not topologically closed after short endpoint collapse";
                    if (!firstWireVertex.IsNull() && !previousEndVertex.IsNull()) {
                        errorMessage += ": "
                            + describeVertexConnection(
                                previousEndVertex,
                                firstWireVertex);
                    }
                    return false;
                }

                BRep_Builder strictWireBuilder;
                TopoDS_Wire strictWire;
                strictWireBuilder.MakeWire(strictWire);
                for (const TopoDS_Edge& rebuiltEdgeUse : rebuiltEdgeUses) {
                    strictWireBuilder.Add(strictWire, rebuiltEdgeUse);
                }
                wireToAdd = strictWire;
                wireToAdd.Orientation(oldWire.Orientation());

                auto validateCanonicalReplacementUses =
                    [&](const TopoDS_Wire& wire, const char* stage) {
                        for (const TopoDS_Edge& expectedUse :
                             expectedReplacementUses) {
                            int sameUseCount = 0;
                            for (TopoDS_Iterator edgeIterator(
                                     wire,
                                     Standard_False,
                                     Standard_False);
                                 edgeIterator.More();
                                 edgeIterator.Next()) {
                                if (edgeIterator.Value().ShapeType() ==
                                        TopAbs_EDGE &&
                                    edgeIterator.Value().IsSame(expectedUse)) {
                                    ++sameUseCount;
                                }
                            }
                            if (sameUseCount != 1) {
                                errorMessage = std::string(stage)
                                    + " does not contain exactly one canonical midline edge use";
                                return false;
                            }
                        }
                        return true;
                    };
                if (!validateCanonicalReplacementUses(
                        wireToAdd,
                        "strict rebuilt wire")) {
                    return false;
                }

                // 局部修复仅处理受影响 wire，最终仍由 ShapeFix_Shape 统一检查和修复整体 shape。
                ShapeFix_Wire wireFix(wireToAdd, forwardFace, faceTolerance);
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
                    errorMessage = "ShapeFix_Wire returned null wire after mapped midline replacement";
                    return false;
                }

                fixedWire.Orientation(oldWire.Orientation());
                wireToAdd = fixedWire;
                if (!validateCanonicalReplacementUses(
                        wireToAdd,
                        "ShapeFix_Wire result")) {
                    return false;
                }
            }
            else {
                const TopoDS_Shape adjustedWireShape = topologyReshaper.Apply(oldWire);
                if (adjustedWireShape.IsNull() || adjustedWireShape.ShapeType() != TopAbs_WIRE) {
                    errorMessage = "failed to propagate midpoint vertices to an unchanged face wire";
                    return false;
                }
                wireToAdd = TopoDS::Wire(adjustedWireShape);
                wireToAdd.Orientation(oldWire.Orientation());
            }

            builder.Add(rebuilt, wireToAdd);
        }

        for (std::size_t taskIndex = 0; taskIndex < replacedTasks.size(); ++taskIndex) {
            if (!replacedTasks[taskIndex]) {
                errorMessage = "face does not contain one of its mapped source edges";
                return false;
            }
        }
        for (std::size_t removalIndex = 0;
             removalIndex < removedTasks.size();
             ++removalIndex) {
            if (!removedTasks[removalIndex]) {
                errorMessage = "face does not contain one of its short endpoint edges";
                return false;
            }
        }

        // 恢复 face 级局部修复，修正 wire 方向以及新中线 edge use 的曲面边界表示。
        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTolerance);
        faceFix.SetMaxTolerance(std::max(faceTolerance, params_.maxGapWidth));
        faceFix.FixWireMode() = 1;
        faceFix.FixOrientationMode() = 1;
        faceFix.FixMissingSeamMode() = 1;
        faceFix.FixAddNaturalBoundMode() = 0;
        faceFix.Perform();

        const TopoDS_Shape fixedShape = faceFix.Result();
        if (!fixedShape.IsNull() && fixedShape.ShapeType() == TopAbs_FACE) {
            rebuiltFace = TopoDS::Face(fixedShape);
        }
        else {
            rebuiltFace = faceFix.Face();
        }
        if (rebuiltFace.IsNull()) {
            errorMessage = "ShapeFix_Face returned null face after mapped midline replacement";
            return false;
        }

        rebuiltFace.Orientation(faceToRebuild.Orientation());
    }
    catch (Standard_Failure& failure) {
        errorMessage = "face rebuild with mapped midline failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFacesByReplacingChain(
    const TopoDS_Shape& currentShape,
    const EdgeChain& keptChain,
    const EdgeChain& replaceChain,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    std::vector<TopoDS_Face> facesToRebuild;
    for (const TopoDS_Face& ownerFace : replaceChain.ownerFaces) {
        if (!ownerFace.IsNull() && !TopologyUtils::IsSameFaceInVector(facesToRebuild, ownerFace)) {
            facesToRebuild.push_back(ownerFace);
        }
    }

    if (facesToRebuild.empty()) {
        errorMessage = "replace chain has no owner faces";
        return false;
    }

    try {
        BRepTools_ReShape reshaper;
        const double tolerance = std::max(params_.maxGapWidth, Precision::Confusion());
        for (const TopoDS_Face& faceToRebuild : facesToRebuild) {
            TopoDS_Face rebuiltFace;
            if (!rebuildFaceByReplacingChainEdges(
                    faceToRebuild,
                    keptChain,
                    replaceChain,
                    tolerance,
                    rebuiltFace,
                    errorMessage)) {
                return false;
            }

            reshaper.Replace(faceToRebuild, rebuiltFace);
        }

        repairedShape = reshaper.Apply(currentShape);
        if (repairedShape.IsNull()) {
            errorMessage = "BRepTools_ReShape produced a null shape";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "shape replacement failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::rebuildFaceByReplacingChainEdges(
    const TopoDS_Face& faceToRebuild,
    const EdgeChain& keptChain,
    const EdgeChain& replaceChain,
    double tolerance,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage) const
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull()) {
        errorMessage = "face to rebuild is null";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "face to rebuild has no support surface";
            return false;
        }

        const double faceTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(faceToRebuild),
            tolerance
        });

        BRep_Builder builder;
        TopoDS_Face rebuilt;
        builder.MakeFace(rebuilt, surface, surfaceLocation, faceTolerance);
        builder.NaturalRestriction(rebuilt, BRep_Tool::NaturalRestriction(forwardFace));
        rebuilt.Orientation(TopAbs_FORWARD);

        bool replacedAnyEdge = false;
        for (TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
             wireExplorer.More();
             wireExplorer.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
            if (oldWire.IsNull()) {
                continue;
            }

            bool wireContainsReplaceEdge = false;
            for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                 edgeExplorer.More();
                 edgeExplorer.Next()) {
                if (findSameEdgeIndex(replaceChain.edges, edgeExplorer.Current()) >= 0) {
                    wireContainsReplaceEdge = true;
                    break;
                }
            }

            TopoDS_Wire wireToAdd = oldWire;
            if (wireContainsReplaceEdge) {
                BRepBuilderAPI_MakeWire wireMaker;
                for (BRepTools_WireExplorer edgeExplorer(oldWire, forwardFace);
                     edgeExplorer.More();
                     edgeExplorer.Next()) {
                    const TopoDS_Edge edgeInWire = edgeExplorer.Current();
                    const int replaceIndex = findSameEdgeIndex(replaceChain.edges, edgeInWire);
                    if (replaceIndex < 0) {
                        wireMaker.Add(edgeInWire);
                        continue;
                    }

                    TopoDS_Edge replacementUse = keptChain.edges[static_cast<std::size_t>(replaceIndex)];
                    replacementUse.Orientation(edgeInWire.Orientation());
                    if (!ensurePCurveOnFace(replacementUse, forwardFace, faceTolerance, errorMessage)) {
                        errorMessage = "failed to repair replacement pcurve: " + errorMessage;
                        return false;
                    }

                    wireMaker.Add(replacementUse);
                    replacedAnyEdge = true;
                }

                if (!wireMaker.IsDone()) {
                    errorMessage = "failed to rebuild wire while replacing chain edges";
                    return false;
                }

                wireToAdd = wireMaker.Wire();
                wireToAdd.Orientation(oldWire.Orientation());

                ShapeFix_Wire wireFix(wireToAdd, forwardFace, faceTolerance);
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
            }

            builder.Add(rebuilt, wireToAdd);
        }

        if (!replacedAnyEdge) {
            errorMessage = "face does not contain replace chain edges";
            return false;
        }

        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTolerance);
        faceFix.SetMaxTolerance(std::max(faceTolerance, params_.maxGapWidth * 10.0));
        faceFix.FixWireMode() = 1;
        faceFix.FixOrientationMode() = 1;
        faceFix.FixMissingSeamMode() = 1;
        faceFix.FixAddNaturalBoundMode() = 0;
        faceFix.Perform();

        const TopoDS_Shape fixedShape = faceFix.Result();
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
        errorMessage = "face rebuild failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

