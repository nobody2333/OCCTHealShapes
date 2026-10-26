#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <TopTools_ListOfShape.hxx>

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

    // 同端点判定只依赖 IsSame（包含 Location，不区分 orientation）。
    // 本轮按无向拓扑端点对索引所有 edge，避免反复扫描不可能匹配的边对。
    // 分组内保持原 edgeFaceMap 索引顺序，外层过滤和日志顺序不变。
    using VertexPair = std::pair<int, int>;
    TopTools_IndexedMapOfShape vertices;
    std::vector<VertexPair> edgeVertexPairs(
        static_cast<std::size_t>(edgeFaceMap.Extent()) + 1);
    std::map<VertexPair, std::vector<int>> edgesByVertexPair;
    for (int edgeIndex = 1; edgeIndex <= edgeFaceMap.Extent(); ++edgeIndex) {
        const TopoDS_Edge edge = TopoDS::Edge(edgeFaceMap.FindKey(edgeIndex));
        TopoDS_Vertex start;
        TopoDS_Vertex end;
        if (!getEdgeVerticesOriented(edge, start, end)) {
            continue;
        }
        const int startIndex = vertices.Add(start);
        const int endIndex = vertices.Add(end);
        const VertexPair pair(
            std::min(startIndex, endIndex), std::max(startIndex, endIndex));
        edgeVertexPairs[static_cast<std::size_t>(edgeIndex)] = pair;
        edgesByVertexPair[pair].push_back(edgeIndex);
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

        const auto matchingEdges = edgesByVertexPair.find(
            edgeVertexPairs[static_cast<std::size_t>(i)]);
        if (matchingEdges == edgesByVertexPair.end()) {
            continue;
        }
        for (const int j : matchingEdges->second) {
            if (j <= i) {
                continue;
            }
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
        // 在双向投影之前提前排除跨 shell 候选。
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

        // 不再要求候选边在公共端点之外相交；通过其余筛选后尝试重建，
        // 是否提交仍由修复后的完整校验决定。

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
                "Crossed-face skipped candidate rejected: "
                + similarityError
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
            minAbsDot = std::min(minAbsDot, GeometryUtils::Clamp01(absDot));
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
        errorMessage = "Curve sampling failed: " + RepairDiagnostics::FailureText(failure);
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
    const double maxGapScore = 1.0 - GeometryUtils::Clamp01(candidate.maxGap / gapDenom);
    const double meanGapScore = 1.0 - GeometryUtils::Clamp01(candidate.meanGap / gapDenom);
    // 长度和切向越接近理想值分数越高。
    const double lengthScore = GeometryUtils::Clamp01(candidate.lengthRatio);
    const double tangentScore = GeometryUtils::Clamp01(candidate.tangentScore);

    // tolA/tolB 用于给容差更稳定的候选略微加权。
    const double tolA = candidate.edgeA.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeA);
    const double tolB = candidate.edgeB.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeB);
    const double tolSpread = std::abs(tolA - tolB) / std::max(std::max(tolA, tolB), Precision::Confusion());
    const double toleranceScore = 1.0 - GeometryUtils::Clamp01(tolSpread * 0.25);

    return 0.35 * maxGapScore +
           0.25 * meanGapScore +
           0.20 * lengthScore +
           0.15 * tangentScore +
           0.05 * toleranceScore;
}

// 选择保留边：明显较小容差优先，否则按参数选择较长边或默认 edgeA。
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
            + RepairDiagnostics::FailureText(failure);
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
            + RepairDiagnostics::FailureText(failure);
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
            + RepairDiagnostics::FailureText(failure);
        return false;
    }
    return std::isfinite(maximumDeviation);
}

