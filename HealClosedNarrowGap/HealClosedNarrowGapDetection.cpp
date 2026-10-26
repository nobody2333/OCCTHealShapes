#include "HealClosedNarrowGap.h"
#include "RepairDiagnostics.h"

#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_POnCurv.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

namespace {
TopoDS_Face firstFaceInList(const TopTools_ListOfShape& shapes)
{
    for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (!face.IsNull()) {
            return face;
        }
    }
    return TopoDS_Face();
}

bool isBoundaryEdgeUsable(const TopoDS_Edge& edge, const TopoDS_Face& ownerFace)
{
    if (edge.IsNull() || ownerFace.IsNull()) {
        return false;
    }

    try {
        if (BRep_Tool::Degenerated(edge) == Standard_True) {
            return false;
        }
        if (BRep_Tool::IsClosed(edge) || BRep_Tool::IsClosed(edge, ownerFace)) {
            return false;
        }
    }
    catch (Standard_Failure&) {
        return false;
    }

    return GeometryUtils::ComputeEdgeLength(edge) > Precision::Confusion();
}

struct BoundaryEdgeUse {
    TopoDS_Edge edge;
    TopoDS_Face ownerFace;
    TopoDS_Vertex firstVertex;
    TopoDS_Vertex lastVertex;
};

struct BoundaryLoop {
    std::vector<TopoDS_Edge> edges;
    std::vector<TopoDS_Face> ownerFaces;
    std::vector<TopoDS_Vertex> vertices;
};

bool edgeHasVertex(const BoundaryEdgeUse& edgeUse, const TopoDS_Vertex& vertex)
{
    return !vertex.IsNull()
        && (edgeUse.firstVertex.IsSame(vertex) || edgeUse.lastVertex.IsSame(vertex));
}

TopoDS_Vertex otherVertexOf(const BoundaryEdgeUse& edgeUse, const TopoDS_Vertex& vertex)
{
    if (edgeUse.firstVertex.IsSame(vertex)) {
        return edgeUse.lastVertex;
    }
    if (edgeUse.lastVertex.IsSame(vertex)) {
        return edgeUse.firstVertex;
    }
    return TopoDS_Vertex();
}

// 仅用于边界分组；集合元素是当前待处理集合中的局部 edge 索引。
struct BoundaryEdgeSets {
    std::vector<std::size_t> parents;
    std::vector<std::size_t> sizes;

    explicit BoundaryEdgeSets(std::size_t count)
        : parents(count), sizes(count, 1)
    {
        std::iota(parents.begin(), parents.end(), std::size_t{0});
    }

    std::size_t find(std::size_t index)
    {
        while (parents[index] != index) {
            parents[index] = parents[parents[index]];
            index = parents[index];
        }
        return index;
    }

    void unite(std::size_t first, std::size_t second)
    {
        first = find(first);
        second = find(second);
        if (first == second) {
            return;
        }
        if (sizes[first] < sizes[second]) {
            std::swap(first, second);
        }
        parents[second] = first;
        sizes[first] += sizes[second];
    }
};

// 输入已由度数分类/并查集保证为单条开放链或简单闭环。
// 排序时直接生成有向 edge use，不在输出后另做一轮拓扑复查。
BoundaryLoop orderBoundaryChain(
    const std::vector<std::size_t>& indices,
    const std::vector<BoundaryEdgeUse>& boundaryEdges,
    const TopTools_IndexedMapOfShape& vertices)
{
    std::map<int, std::vector<std::size_t>> incident;
    for (const std::size_t index : indices) {
        const BoundaryEdgeUse& edge = boundaryEdges[index];
        incident[vertices.FindIndex(edge.firstVertex)].push_back(index);
        incident[vertices.FindIndex(edge.lastVertex)].push_back(index);
    }

    int start = vertices.FindIndex(boundaryEdges[indices.front()].firstVertex);
    for (const auto& entry : incident) {
        assert(entry.second.size() == 1 || entry.second.size() == 2);
        if (entry.second.size() == 1) {
            start = entry.first;
            break;
        }
    }

    BoundaryLoop chain;
    TopoDS_Vertex current = TopoDS::Vertex(vertices.FindKey(start));
    chain.vertices.push_back(current);
    std::size_t previous = boundaryEdges.size();
    while (true) {
        const auto& adjacent = incident.at(vertices.FindIndex(current));
        std::size_t next = adjacent.front();
        if (next == previous) {
            if (adjacent.size() == 1) {
                break;
            }
            next = adjacent.back();
        }

        const BoundaryEdgeUse& edge = boundaryEdges[next];
        assert(edgeHasVertex(edge, current));
        const TopoDS_Vertex finish = otherVertexOf(edge, current);
        chain.edges.push_back(TopologyUtils::OrientEdgeFromTo(edge.edge, current, finish));
        chain.ownerFaces.push_back(edge.ownerFace);
        chain.vertices.push_back(finish);
        previous = next;
        current = finish;
        if (current.IsSame(chain.vertices.front())) {
            break;
        }
    }
    assert(chain.edges.size() == indices.size());
    assert(chain.vertices.size() == chain.edges.size() + 1);
    return chain;
}

std::vector<BoundaryLoop> extractBoundaryLoops(
    const std::vector<BoundaryEdgeUse>& boundaryEdges)
{
    std::vector<BoundaryLoop> loops;
    if (boundaryEdges.empty()) {
        return loops;
    }

    // IndexedMap 使用拓扑身份（包括 location，忽略 orientation），不合并近邻坐标。
    TopTools_IndexedMapOfShape vertices;
    for (const BoundaryEdgeUse& edge : boundaryEdges) {
        vertices.Add(edge.firstVertex);
        vertices.Add(edge.lastVertex);
    }

    std::vector<std::size_t> allIndices(boundaryEdges.size());
    std::iota(allIndices.begin(), allIndices.end(), std::size_t{0});
    std::vector<std::vector<std::size_t>> pending;
    pending.push_back(std::move(allIndices));
    std::set<std::vector<std::size_t>> emitted;
    auto emitLoop = [&](std::vector<std::size_t> indices, BoundaryLoop loop) {
        if (indices.size() < 3) {
            return false;
        }
        std::sort(indices.begin(), indices.end());
        if (emitted.insert(std::move(indices)).second) {
            loops.push_back(std::move(loop));
        }
        return true;
    };

    while (!pending.empty()) {
        std::vector<std::size_t> indices = std::move(pending.back());
        pending.pop_back();

        std::map<int, std::vector<std::size_t>> incident;
        for (std::size_t local = 0; local < indices.size(); ++local) {
            const BoundaryEdgeUse& edge = boundaryEdges[indices[local]];
            incident[vertices.FindIndex(edge.firstVertex)].push_back(local);
            incident[vertices.FindIndex(edge.lastVertex)].push_back(local);
        }

        auto groupEdges = [&](bool joinBranches) {
            BoundaryEdgeSets sets(indices.size());
            for (const auto& entry : incident) {
                const auto& adjacent = entry.second;
                if (joinBranches || adjacent.size() == 2) {
                    for (std::size_t i = 1; i < adjacent.size(); ++i) {
                        sets.unite(adjacent.front(), adjacent[i]);
                    }
                }
            }
            std::map<std::size_t, std::vector<std::size_t>> groups;
            for (std::size_t local = 0; local < indices.size(); ++local) {
                groups[sets.find(local)].push_back(indices[local]);
            }
            // 按原始最小边索引保持分量/链的确定性顺序。
            std::vector<std::vector<std::size_t>> result;
            for (auto& entry : groups) {
                result.push_back(std::move(entry.second));
            }
            std::sort(result.begin(), result.end(),
                [](const auto& first, const auto& second) {
                    return first.front() < second.front();
                });
            return result;
        };

        // 必须在任何顺序遍历前分离连通分量；剩余边也重新走此步骤。
        auto components = groupEdges(true);
        if (components.size() > 1) {
            for (auto it = components.rbegin(); it != components.rend(); ++it) {
                pending.push_back(std::move(*it));
            }
            continue;
        }

        bool branched = false;
        bool open = false;
        for (const auto& entry : incident) {
            branched = branched || entry.second.size() > 2;
            open = open || entry.second.size() == 1;
        }
        if (!branched) {
            if (!open) {
                emitLoop(indices, orderBoundaryChain(indices, boundaryEdges, vertices));
            }
            continue;
        }

        // 只在度为 2 的顶点合并边：每组必为链，或回到同一分叉点的简单环。
        const auto groups = groupEdges(false);
        std::vector<BoundaryLoop> chains;
        std::map<std::pair<int, int>, std::vector<std::size_t>> endpointPairs;
        std::set<std::size_t> covered;
        for (std::size_t group = 0; group < groups.size(); ++group) {
            chains.push_back(orderBoundaryChain(groups[group], boundaryEdges, vertices));
            const BoundaryLoop& chain = chains.back();
            int first = vertices.FindIndex(chain.vertices.front());
            int last = vertices.FindIndex(chain.vertices.back());
            if (first == last) {
                if (emitLoop(groups[group], chain)) {
                    covered.insert(groups[group].begin(), groups[group].end());
                }
            }
            else {
                if (first > last) {
                    std::swap(first, last);
                }
                endpointPairs[{first, last}].push_back(group);
            }
        }

        // 本轮完整配对，不因先前组合覆盖了某条链而跳过其他链对。
        for (const auto& entry : endpointPairs) {
            const auto& matching = entry.second;
            for (std::size_t i = 0; i < matching.size(); ++i) {
                for (std::size_t j = i + 1; j < matching.size(); ++j) {
                    const std::size_t a = matching[i], b = matching[j];
                    if (groups[a].size() + groups[b].size() < 3) {
                        continue;
                    }
                    BoundaryLoop loop = chains[a];
                    const BoundaryLoop& second = chains[b];
                    // 开放链均从拓扑索引较小的端点开始，第二条需反向追加。
                    assert(loop.vertices.front().IsSame(second.vertices.front()));
                    for (std::size_t k = second.edges.size(); k > 0; --k) {
                        loop.edges.push_back(TopologyUtils::OrientEdgeFromTo(
                            second.edges[k - 1], second.vertices[k], second.vertices[k - 1]));
                        loop.ownerFaces.push_back(second.ownerFaces[k - 1]);
                        loop.vertices.push_back(second.vertices[k - 1]);
                    }
                    std::vector<std::size_t> loopIndices = groups[a];
                    loopIndices.insert(loopIndices.end(), groups[b].begin(), groups[b].end());
                    if (emitLoop(loopIndices, std::move(loop))) {
                        covered.insert(loopIndices.begin(), loopIndices.end());
                    }
                }
            }
        }

        // 分组由当前边集合唯一决定：没有移除边，下一轮分组也不会改变。
        if (covered.empty()) {
            continue;
        }
        std::vector<std::size_t> remaining;
        for (const std::size_t index : indices) {
            if (covered.find(index) == covered.end()) {
                remaining.push_back(index);
            }
        }
        if (!remaining.empty()) {
            pending.push_back(std::move(remaining));
        }
    }
    return loops;
}

std::vector<BoundaryLoop> buildBoundaryLoops(const TopoDS_Shape& shape)
{
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);

    std::vector<BoundaryEdgeUse> boundaryEdges;
    for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
        if (edgeFaceMap.FindFromIndex(i).Extent() != 1) {
            continue;
        }

        const TopoDS_Edge edge = TopoDS::Edge(edgeFaceMap.FindKey(i));
        const TopoDS_Face ownerFace = firstFaceInList(edgeFaceMap.FindFromIndex(i));
        if (!isBoundaryEdgeUsable(edge, ownerFace)) {
            continue;
        }

        BoundaryEdgeUse edgeUse;
        edgeUse.edge = edge;
        edgeUse.ownerFace = ownerFace;
        if (!TopologyUtils::GetEdgeVertices(edge, edgeUse.firstVertex, edgeUse.lastVertex)) {
            continue;
        }
        boundaryEdges.push_back(edgeUse);
    }

    return extractBoundaryLoops(boundaryEdges);
}
} // namespace

std::vector<HealClosedNarrowGap::ClosedGapCandidate> HealClosedNarrowGap::extractCandidates(
    const TopoDS_Shape& shape) const
{
    std::vector<ClosedGapCandidate> candidates;
    if (shape.IsNull()) {
        return candidates;
    }

    const std::vector<BoundaryLoop> loops = buildBoundaryLoops(shape);
    for (const BoundaryLoop& loop : loops) {
        const int edgeCount = static_cast<int>(loop.edges.size());
        if (edgeCount < 3) {
            continue;
        }

        std::vector<double> edgeLengths(static_cast<std::size_t>(edgeCount), 0.0);
        bool edgeLengthsValid = true;
        for (int edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
            const double length = GeometryUtils::ComputeEdgeLength(
                loop.edges[static_cast<std::size_t>(edgeIndex)]);
            if (!std::isfinite(length) || length <= Precision::Confusion()) {
                edgeLengthsValid = false;
                break;
            }
            edgeLengths[static_cast<std::size_t>(edgeIndex)] = length;
        }
        if (!edgeLengthsValid) {
            continue;
        }

        std::vector<int> tipIndices;
        std::vector<double> tipScores(edgeCount, 0.0);
        for (int i = 0; i < edgeCount; ++i) {
            const int prevIndex = (i + edgeCount - 1) % edgeCount;
            const double angle = GeometryUtils::ComputeIntersectionAngle(
                loop.edges[prevIndex],
                loop.edges[i],
                loop.vertices[static_cast<std::size_t>(i)]);
            if (angle < 0.0 || angle > params_.maxTipAngle) {
                continue;
            }

            tipIndices.push_back(i);
            tipScores[static_cast<std::size_t>(i)] =
                1.0 - GeometryUtils::Clamp01(angle / std::max(params_.maxTipAngle, Precision::Confusion()));
        }

        auto buildForwardChain = [&](int startVertexIndex,
                                     int endVertexIndex,
                                     EdgeChain& chain) {
            chain = EdgeChain();
            startVertexIndex = (startVertexIndex + edgeCount) % edgeCount;
            endVertexIndex = (endVertexIndex + edgeCount) % edgeCount;
            chain.vertices.push_back(
                loop.vertices[static_cast<std::size_t>(startVertexIndex)]);

            int edgeIndex = startVertexIndex;
            int traversed = 0;
            while (edgeIndex != endVertexIndex && traversed < edgeCount) {
                const int nextVertexIndex = (edgeIndex + 1) % edgeCount;
                const TopoDS_Vertex from =
                    loop.vertices[static_cast<std::size_t>(edgeIndex)];
                const TopoDS_Vertex to =
                    loop.vertices[static_cast<std::size_t>(nextVertexIndex)];
                const TopoDS_Edge edge = TopologyUtils::OrientEdgeFromTo(
                    loop.edges[static_cast<std::size_t>(edgeIndex)],
                    from,
                    to);
                chain.edges.push_back(edge);
                chain.ownerFaces.push_back(
                    loop.ownerFaces[static_cast<std::size_t>(edgeIndex)]);
                chain.vertices.push_back(to);
                chain.length +=
                    edgeLengths[static_cast<std::size_t>(edgeIndex)];
                edgeIndex = nextVertexIndex;
                ++traversed;
            }

            return edgeIndex == endVertexIndex &&
                !chain.edges.empty() &&
                chain.vertices.size() == chain.edges.size() + 1 &&
                std::isfinite(chain.length) &&
                chain.length > Precision::Confusion();
        };

        auto reverseChain = [&](const EdgeChain& source, EdgeChain& reversed) {
            reversed = EdgeChain();
            if (source.edges.empty() ||
                source.vertices.size() != source.edges.size() + 1) {
                return false;
            }

            reversed.vertices.push_back(source.vertices.back());
            for (std::size_t reverseIndex = source.edges.size();
                 reverseIndex > 0;
                 --reverseIndex) {
                const std::size_t sourceIndex = reverseIndex - 1;
                const TopoDS_Vertex from = source.vertices[sourceIndex + 1];
                const TopoDS_Vertex to = source.vertices[sourceIndex];
                reversed.edges.push_back(TopologyUtils::OrientEdgeFromTo(
                    source.edges[sourceIndex], from, to));
                reversed.ownerFaces.push_back(
                    source.ownerFaces[sourceIndex]);
                reversed.vertices.push_back(to);
            }
            reversed.length = source.length;
            return true;
        };

        auto makeTipEnd = [&](int vertexIndex) {
            GapEnd end;
            end.type = GapEndType::SharedTip;
            end.vertexA = loop.vertices[static_cast<std::size_t>(vertexIndex)];
            end.vertexB = end.vertexA;
            end.score = tipScores[static_cast<std::size_t>(vertexIndex)];
            return end;
        };

        auto makeShortEdgeEnd = [&](int edgeIndex,
                                    const TopoDS_Vertex& vertexA,
                                    const TopoDS_Vertex& vertexB,
                                    double score) {
            GapEnd end;
            end.type = GapEndType::ShortEdge;
            end.vertexA = vertexA;
            end.vertexB = vertexB;
            end.connectorEdge =
                loop.edges[static_cast<std::size_t>(edgeIndex)];
            end.connectorOwnerFace =
                loop.ownerFaces[static_cast<std::size_t>(edgeIndex)];
            end.score = score;
            return end;
        };

        auto appendCandidate = [&](const EdgeChain& chainA,
                                   const EdgeChain& chainB,
                                   const GapEnd& startEnd,
                                   const GapEnd& finishEnd) {
            if (chainA.edges.empty() || chainB.edges.empty() ||
                chainA.vertices.size() != chainA.edges.size() + 1 ||
                chainB.vertices.size() != chainB.edges.size() + 1 ||
                !isCandidateEndConsistent(
                    startEnd,
                    chainA.vertices.front(),
                    chainB.vertices.front()) ||
                !isCandidateEndConsistent(
                    finishEnd,
                    chainA.vertices.back(),
                    chainB.vertices.back())) {
                return;
            }

            ClosedGapCandidate candidate;
            candidate.chainA = chainA;
            candidate.chainB = chainB;
            candidate.startEnd = startEnd;
            candidate.finishEnd = finishEnd;
            candidate.endpointScore =
                0.5 * (startEnd.score + finishEnd.score);

            // 长度比只依赖已计算的 chain 长度，先做廉价过滤，避免无效组合进入双向投影。
            const double longerChainLength = std::max(chainA.length, chainB.length);
            if (!std::isfinite(longerChainLength) ||
                longerChainLength <= Precision::Confusion()) {
                return;
            }
            candidate.lengthRatio =
                std::min(chainA.length, chainB.length) / longerChainLength;
            if (!std::isfinite(candidate.lengthRatio) ||
                candidate.lengthRatio + Precision::Confusion() <
                    params_.minLengthRatio) {
                return;
            }

            auto endpointPoint = [](const GapEnd& end, gp_Pnt& point) {
                if (end.vertexA.IsNull() || end.vertexB.IsNull()) {
                    return false;
                }
                try {
                    const gp_Pnt pointA = BRep_Tool::Pnt(end.vertexA);
                    const gp_Pnt pointB = BRep_Tool::Pnt(end.vertexB);
                    point.SetCoord(
                        0.5 * (pointA.X() + pointB.X()),
                        0.5 * (pointA.Y() + pointB.Y()),
                        0.5 * (pointA.Z() + pointB.Z()));
                    return GeometryUtils::IsFinitePoint(point);
                }
                catch (const Standard_Failure&) {
                    return false;
                }
            };

            gp_Pnt startPoint;
            gp_Pnt finishPoint;
            if (!endpointPoint(startEnd, startPoint) ||
                !endpointPoint(finishEnd, finishPoint)) {
                return;
            }
            candidate.endpointSpan = startPoint.Distance(finishPoint);
            if (!std::isfinite(candidate.endpointSpan) ||
                candidate.endpointSpan <= Precision::Confusion()) {
                return;
            }

            std::string gapError;
            if (!evaluateBidirectionalGap(
                    candidate.chainA,
                    candidate.chainB,
                    candidate.maxWidthAToB,
                    candidate.maxWidthBToA,
                    candidate.meanWidth,
                    candidate.usedWidthSampleCount,
                    gapError)) {
                return;
            }

            candidate.maxWidth = std::max(
                candidate.maxWidthAToB,
                candidate.maxWidthBToA);
            candidate.lengthWidthRatio =
                0.5 * (chainA.length + chainB.length) /
                std::max(candidate.meanWidth, Precision::Confusion());

            double endpointWidthReference = candidate.maxWidth;
            try {
                if (startEnd.type == GapEndType::ShortEdge) {
                    endpointWidthReference = std::max(
                        endpointWidthReference,
                        BRep_Tool::Pnt(startEnd.vertexA).Distance(
                            BRep_Tool::Pnt(startEnd.vertexB)));
                }
                if (finishEnd.type == GapEndType::ShortEdge) {
                    endpointWidthReference = std::max(
                        endpointWidthReference,
                        BRep_Tool::Pnt(finishEnd.vertexA).Distance(
                            BRep_Tool::Pnt(finishEnd.vertexB)));
                }
            }
            catch (const Standard_Failure&) {
                return;
            }
            candidate.endpointSpanWidthRatio =
                candidate.endpointSpan /
                std::max(endpointWidthReference, Precision::Confusion());
            const bool hasShortEdgeEnd =
                startEnd.type == GapEndType::ShortEdge ||
                finishEnd.type == GapEndType::ShortEdge;

            if (candidate.maxWidth > params_.maxGapWidth ||
                candidate.lengthWidthRatio + Precision::Confusion() <
                    params_.minLengthWidthRatio ||
                (hasShortEdgeEnd &&
                    candidate.endpointSpanWidthRatio +
                            Precision::Confusion() <
                        params_.minLengthWidthRatio)) {
                return;
            }

            const double widthScore =
                1.0 - GeometryUtils::Clamp01(candidate.maxWidth / params_.maxGapWidth);
            const double ratioScore = GeometryUtils::Clamp01(candidate.lengthRatio);
            const double slenderScore = GeometryUtils::Clamp01(
                candidate.lengthWidthRatio /
                std::max(params_.minLengthWidthRatio * 2.0, 1.0));
            candidate.confidence =
                0.40 * widthScore +
                0.25 * slenderScore +
                0.20 * ratioScore +
                0.15 * GeometryUtils::Clamp01(candidate.endpointScore);
            candidates.push_back(candidate);
        };

        // 原有双尖点路径保持不变：两个尖点把闭合 loop 分成方向一致的两侧 chain。
        if (tipIndices.size() >= 2) {
            for (std::size_t lhs = 0; lhs < tipIndices.size(); ++lhs) {
                for (std::size_t rhs = lhs + 1; rhs < tipIndices.size(); ++rhs) {
                    const int startIndex = tipIndices[lhs];
                    const int endIndex = tipIndices[rhs];
                    EdgeChain chainA;
                    EdgeChain reversePath;
                    EdgeChain chainB;
                    if (!buildForwardChain(startIndex, endIndex, chainA) ||
                        !buildForwardChain(endIndex, startIndex, reversePath) ||
                        !reverseChain(reversePath, chainB)) {
                        continue;
                    }
                    appendCandidate(
                        chainA,
                        chainB,
                        makeTipEnd(startIndex),
                        makeTipEnd(endIndex));
                }
            }
            continue;
        }

        if (!params_.detectShortEdgeEnds) {
            continue;
        }

        // 非尖点端部必须是有效短 edge，且它两侧相邻 edge 都明显更长。
        const double maxEndEdgeLength =
            params_.maxGapWidth * params_.maxEndEdgeLengthFactor;
        std::vector<int> shortEndEdgeIndices;
        std::vector<double> shortEndScores(
            static_cast<std::size_t>(edgeCount), 0.0);
        for (int edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
            const int previousIndex =
                (edgeIndex + edgeCount - 1) % edgeCount;
            const int nextIndex = (edgeIndex + 1) % edgeCount;
            const double edgeLength =
                edgeLengths[static_cast<std::size_t>(edgeIndex)];
            const double shorterAdjacentLength = std::min(
                edgeLengths[static_cast<std::size_t>(previousIndex)],
                edgeLengths[static_cast<std::size_t>(nextIndex)]);
            if (edgeLength > maxEndEdgeLength ||
                shorterAdjacentLength + Precision::Confusion() <
                    edgeLength * params_.minAdjacentEdgeLengthRatio) {
                continue;
            }

            shortEndEdgeIndices.push_back(edgeIndex);
            shortEndScores[static_cast<std::size_t>(edgeIndex)] =
                1.0 - GeometryUtils::Clamp01(
                    edgeLength /
                    std::max(maxEndEdgeLength, Precision::Confusion()));
        }

        if (tipIndices.size() == 1) {
            // 一个尖点和一条对侧短端边共同定义两侧 chain。
            const int tipIndex = tipIndices.front();
            for (const int endEdgeIndex : shortEndEdgeIndices) {
                EdgeChain chainA;
                EdgeChain reversePath;
                EdgeChain chainB;
                if (!buildForwardChain(tipIndex, endEdgeIndex, chainA) ||
                    !buildForwardChain(
                        (endEdgeIndex + 1) % edgeCount,
                        tipIndex,
                        reversePath) ||
                    !reverseChain(reversePath, chainB)) {
                    continue;
                }

                const GapEnd shortEnd = makeShortEdgeEnd(
                    endEdgeIndex,
                    chainA.vertices.back(),
                    chainB.vertices.back(),
                    shortEndScores[static_cast<std::size_t>(endEdgeIndex)]);
                appendCandidate(
                    chainA,
                    chainB,
                    makeTipEnd(tipIndex),
                    shortEnd);
            }
            continue;
        }

        // 没有尖点时，两条短端边分别收缩为中线的起点和终点。
        for (std::size_t lhs = 0; lhs < shortEndEdgeIndices.size(); ++lhs) {
            for (std::size_t rhs = lhs + 1;
                 rhs < shortEndEdgeIndices.size();
                 ++rhs) {
                const int startEdgeIndex = shortEndEdgeIndices[lhs];
                const int finishEdgeIndex = shortEndEdgeIndices[rhs];
                EdgeChain chainA;
                EdgeChain reversePath;
                EdgeChain chainB;
                if (!buildForwardChain(
                        (startEdgeIndex + 1) % edgeCount,
                        finishEdgeIndex,
                        chainA) ||
                    !buildForwardChain(
                        (finishEdgeIndex + 1) % edgeCount,
                        startEdgeIndex,
                        reversePath) ||
                    !reverseChain(reversePath, chainB)) {
                    continue;
                }

                const GapEnd startEnd = makeShortEdgeEnd(
                    startEdgeIndex,
                    chainA.vertices.front(),
                    chainB.vertices.front(),
                    shortEndScores[static_cast<std::size_t>(startEdgeIndex)]);
                const GapEnd finishEnd = makeShortEdgeEnd(
                    finishEdgeIndex,
                    chainA.vertices.back(),
                    chainB.vertices.back(),
                    shortEndScores[static_cast<std::size_t>(finishEdgeIndex)]);
                appendCandidate(chainA, chainB, startEnd, finishEnd);
            }
        }
    }

    return candidates;
}

bool HealClosedNarrowGap::sampleChainPoint(const EdgeChain& chain, double normalized, gp_Pnt& point) const
{
    if (chain.edges.empty() || chain.vertices.size() != chain.edges.size() + 1 ||
        !std::isfinite(chain.length) ||
        chain.length <= Precision::Confusion()) {
        return false;
    }

    try {
        const double targetLength = chain.length * GeometryUtils::Clamp01(normalized);
        double accumulated = 0.0;
        for (std::size_t i = 0; i < chain.edges.size(); ++i) {
            const TopoDS_Edge& edge = chain.edges[i];
            const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
            if (!std::isfinite(edgeLength) ||
                edgeLength <= Precision::Confusion()) {
                return false;
            }

            const bool isLastEdge = i + 1 == chain.edges.size();
            if (!isLastEdge &&
                accumulated + edgeLength < targetLength) {
                accumulated += edgeLength;
                continue;
            }

            const double localLength = std::max(
                0.0,
                std::min(edgeLength, targetLength - accumulated));
            const double startParameter =
                BRep_Tool::Parameter(chain.vertices[i], edge);
            const double endParameter =
                BRep_Tool::Parameter(chain.vertices[i + 1], edge);
            if (!std::isfinite(startParameter) ||
                !std::isfinite(endParameter) ||
                std::abs(endParameter - startParameter) <=
                    Precision::PConfusion()) {
                return false;
            }

            // BRepAdaptor_Curve 会应用 edge 的 TopLoc_Location，投影点与源点位于同一全局坐标系。
            BRepAdaptor_Curve curve(edge);
            double parameter = startParameter;
            if (localLength >= edgeLength - Precision::Confusion()) {
                parameter = endParameter;
            }
            else if (localLength > Precision::Confusion()) {
                const double initialParameter =
                    startParameter
                    + (endParameter - startParameter)
                        * (localLength / edgeLength);
                const double signedAbscissa =
                    startParameter <= endParameter
                        ? localLength
                        : -localLength;
                GCPnts_AbscissaPoint abscissaPoint(
                    curve,
                    signedAbscissa,
                    startParameter,
                    initialParameter,
                    Precision::Confusion());
                if (!abscissaPoint.IsDone()) {
                    return false;
                }
                parameter = abscissaPoint.Parameter();
            }

            point = curve.Value(parameter);
            return GeometryUtils::IsFinitePoint(point);
        }
    }
    catch (const Standard_Failure&) {
        return false;
    }

    return false;
}

bool HealClosedNarrowGap::findNearestPointOnChain(
    const gp_Pnt& sourcePoint,
    const EdgeChain& targetChain,
    std::size_t& targetEdgeIndex,
    double& targetParameter,
    gp_Pnt& targetPoint,
    double& distance,
    std::string& errorMessage) const
{
    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    targetEdgeIndex = invalidIndex;
    targetParameter = 0.0;
    targetPoint = gp_Pnt();
    distance = std::numeric_limits<double>::infinity();

    if (!GeometryUtils::IsFinitePoint(sourcePoint) ||
        targetChain.edges.empty() ||
        targetChain.vertices.size() != targetChain.edges.size() + 1) {
        errorMessage = "nearest chain projection failure: target chain data is invalid";
        return false;
    }

    try {
        for (std::size_t edgeIndex = 0;
             edgeIndex < targetChain.edges.size();
             ++edgeIndex) {
            const TopoDS_Edge& edge = targetChain.edges[edgeIndex];
            if (edge.IsNull() ||
                targetChain.vertices[edgeIndex].IsNull() ||
                targetChain.vertices[edgeIndex + 1].IsNull()) {
                errorMessage = "nearest chain projection failure: target chain contains null topology";
                return false;
            }

            BRepAdaptor_Curve curve(edge);
            const double startParameter =
                BRep_Tool::Parameter(targetChain.vertices[edgeIndex], edge);
            const double endParameter =
                BRep_Tool::Parameter(targetChain.vertices[edgeIndex + 1], edge);
            const double lowerParameter = std::max(
                curve.FirstParameter(),
                std::min(startParameter, endParameter));
            const double upperParameter = std::min(
                curve.LastParameter(),
                std::max(startParameter, endParameter));
            if (!std::isfinite(startParameter) ||
                !std::isfinite(endParameter) ||
                !std::isfinite(lowerParameter) ||
                !std::isfinite(upperParameter) ||
                upperParameter - lowerParameter <= Precision::PConfusion()) {
                errorMessage = "nearest chain projection failure: target edge has an invalid trimmed range";
                return false;
            }

            auto considerProjection =
                [&](double parameter, const gp_Pnt& candidatePoint) {
                    if (!std::isfinite(parameter) ||
                        !GeometryUtils::IsFinitePoint(candidatePoint) ||
                        parameter < lowerParameter - Precision::PConfusion() ||
                        parameter > upperParameter + Precision::PConfusion()) {
                        return;
                    }

                    const double candidateDistance =
                        sourcePoint.Distance(candidatePoint);
                    if (!std::isfinite(candidateDistance)) {
                        return;
                    }
                    if (targetEdgeIndex == invalidIndex ||
                        candidateDistance + Precision::Confusion() < distance) {
                        distance = candidateDistance;
                        targetEdgeIndex = edgeIndex;
                        targetParameter = parameter;
                        targetPoint = candidatePoint;
                    }
                };

            // 极值求解限制在实际裁剪参数范围内，并显式加入两个端点作为兜底候选。
            considerProjection(lowerParameter, curve.Value(lowerParameter));
            considerProjection(upperParameter, curve.Value(upperParameter));
            Extrema_ExtPC extrema(
                sourcePoint,
                curve,
                lowerParameter,
                upperParameter,
                Precision::Confusion());
            if (extrema.IsDone()) {
                for (int solutionIndex = 1;
                     solutionIndex <= extrema.NbExt();
                     ++solutionIndex) {
                    const Extrema_POnCurv& pointOnCurve =
                        extrema.Point(solutionIndex);
                    considerProjection(
                        pointOnCurve.Parameter(),
                        pointOnCurve.Value());
                }
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "nearest chain projection failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    if (targetEdgeIndex == invalidIndex ||
        !std::isfinite(targetParameter) ||
        !GeometryUtils::IsFinitePoint(targetPoint) ||
        !std::isfinite(distance)) {
        errorMessage = "nearest chain projection failure: no valid projection on trimmed target edges";
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::evaluateBidirectionalGap(
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    double& maxWidthAToB,
    double& maxWidthBToA,
    double& meanWidth,
    int& usedSampleCount,
    std::string& errorMessage) const
{
    maxWidthAToB = 0.0;
    maxWidthBToA = 0.0;
    meanWidth = 0.0;
    usedSampleCount = computeWidthSampleCount(
        chainA.length,
        chainB.length);
    const int sampleCount = usedSampleCount;
    double distanceSum = 0.0;
    int distanceCount = 0;

    auto accumulateDirectedDistances =
        [&](const EdgeChain& sourceChain,
            const EdgeChain& targetChain,
            double& directedMaximum,
            const char* directionName) {
            for (int sampleIndex = 0;
                 sampleIndex < sampleCount;
                 ++sampleIndex) {
                const double normalized =
                    static_cast<double>(sampleIndex)
                    / static_cast<double>(sampleCount - 1);
                gp_Pnt sourcePoint;
                if (!sampleChainPoint(
                        sourceChain,
                        normalized,
                        sourcePoint)) {
                    errorMessage = std::string(
                        "bidirectional gap evaluation failure: failed to sample ")
                        + directionName;
                    return false;
                }

                std::size_t targetEdgeIndex = 0;
                double targetParameter = 0.0;
                gp_Pnt targetPoint;
                double nearestDistance = 0.0;
                std::string projectionError;
                if (!findNearestPointOnChain(
                        sourcePoint,
                        targetChain,
                        targetEdgeIndex,
                        targetParameter,
                        targetPoint,
                        nearestDistance,
                        projectionError)) {
                    errorMessage = "bidirectional gap evaluation failure: "
                        + projectionError;
                    return false;
                }

                // 最大宽度取所有双向采样距离的最大值；任一点超限后，后续采样不可能使候选重新合格。
                if (nearestDistance > params_.maxGapWidth) {
                    errorMessage = std::string(
                        "bidirectional gap evaluation rejected: ")
                        + directionName
                        + " sample distance exceeds maxGapWidth";
                    return false;
                }

                directedMaximum = std::max(
                    directedMaximum,
                    nearestDistance);
                distanceSum += nearestDistance;
                ++distanceCount;
            }
            return true;
        };

    // 两个方向都必须覆盖，避免只检查 A 靠近 B，却遗漏 B 的局部外张。
    if (!accumulateDirectedDistances(
            chainA,
            chainB,
            maxWidthAToB,
            "chainA") ||
        !accumulateDirectedDistances(
            chainB,
            chainA,
            maxWidthBToA,
            "chainB")) {
        return false;
    }

    if (distanceCount != 2 * sampleCount ||
        !std::isfinite(distanceSum) ||
        !std::isfinite(maxWidthAToB) ||
        !std::isfinite(maxWidthBToA)) {
        errorMessage = "bidirectional gap evaluation failure: distance statistics are invalid";
        return false;
    }

    meanWidth = distanceSum / static_cast<double>(distanceCount);
    if (!std::isfinite(meanWidth)) {
        errorMessage = "bidirectional gap evaluation failure: mean distance is invalid";
        return false;
    }
    return true;
}

