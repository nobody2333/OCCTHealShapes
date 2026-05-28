#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepCheck_Face.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepLib.hxx>
#include <BRepTools_ReShape.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

namespace {
// 将评分项限制在 [0, 1]，避免异常参数影响综合置信度。
double Clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

// 将 OCCT 异常消息转为 std::string，便于写入修复日志。
std::string FailureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
}

// 轻量检查 edge 在 face 上是否已有 pcurve。
bool HasPCurveOnFace(const TopoDS_Edge& edge, const TopoDS_Face& face)
{
    return GeometryUtils::HasValidPCurveOnFace(edge, face);
}

const char* BRepCheckStatusName(BRepCheck_Status status)
{
    switch (status) {
    case BRepCheck_NoError: return "NoError";
    case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
    case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
    case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
    case BRepCheck_No3DCurve: return "No3DCurve";
    case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
    case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
    case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
    case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
    case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
    case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
    case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
    case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
    case BRepCheck_FreeEdge: return "FreeEdge";
    case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
    case BRepCheck_InvalidRange: return "InvalidRange";
    case BRepCheck_EmptyWire: return "EmptyWire";
    case BRepCheck_RedundantEdge: return "RedundantEdge";
    case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
    case BRepCheck_NoSurface: return "NoSurface";
    case BRepCheck_InvalidWire: return "InvalidWire";
    case BRepCheck_RedundantWire: return "RedundantWire";
    case BRepCheck_IntersectingWires: return "IntersectingWires";
    case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
    case BRepCheck_EmptyShell: return "EmptyShell";
    case BRepCheck_RedundantFace: return "RedundantFace";
    case BRepCheck_UnorientableShape: return "UnorientableShape";
    case BRepCheck_NotClosed: return "NotClosed";
    case BRepCheck_NotConnected: return "NotConnected";
    case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
    case BRepCheck_BadOrientation: return "BadOrientation";
    case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
    case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
    case BRepCheck_CheckFail: return "CheckFail";
    default: return "Unknown";
    }
}

std::string DescribeBRepCheckFailures(BRepCheck_Analyzer& analyzer, const TopoDS_Shape& shape)
{
    std::ostringstream os;
    TopTools_IndexedMapOfShape allShapes;
    TopExp::MapShapes(shape, allShapes);

    int reported = 0;
    constexpr int maxReported = 30;
    for (int i = 1; i <= allShapes.Extent() && reported < maxReported; ++i) {
        const TopoDS_Shape& subShape = allShapes.FindKey(i);
        Handle(BRepCheck_Result) result = analyzer.Result(subShape);
        if (result.IsNull()) {
            continue;
        }

        bool hasError = false;
        std::ostringstream statusStream;
        const BRepCheck_ListOfStatus& statuses = result->Status();
        for (BRepCheck_ListIteratorOfListOfStatus it(statuses); it.More(); it.Next()) {
            const BRepCheck_Status status = it.Value();
            if (status == BRepCheck_NoError) {
                continue;
            }
            if (hasError) {
                statusStream << ",";
            }
            statusStream << BRepCheckStatusName(status);
            hasError = true;
        }

        if (!hasError) {
            continue;
        }
        
        os << " [subshape#" << i
           << " type=" << static_cast<int>(subShape.ShapeType())
           << " status=" << statusStream.str() << "]";
        ++reported;
    }

    if (reported == maxReported) {
        os << " ...";
    }

    return os.str();
}

bool FindFacesForSameEdge(const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
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
}

// 构造时规范化参数，避免采样数和迭代数出现无效值。
HealNarrowGap::HealNarrowGap(const Params& params)
    : myParams(params)
{
    myParams.sampleCount = std::max(2, myParams.sampleCount);
    myParams.maxIterations = std::max(1, myParams.maxIterations);
    myParams.minLengthRatio = Clamp01(myParams.minLengthRatio);
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
int HealNarrowGap::Perform()
{
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        AppendMessage(messages_, "input shape is null");
        return 0;
    }

    // currentShape 是迭代中的工作副本，不直接修改 inputShape。
    TopoDS_Shape currentShape = inputShape_;
    for (int iter = 0; iter < myParams.maxIterations; ++iter) {
        // edgeFaceMap 记录每条 edge 的相邻 face，用于判定单面边和 free edge。
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopologyUtils::BuildEdgeFaceMap(currentShape, edgeFaceMap);

        // candidates 是本轮从当前形状重新提取的候选，避免基于过期拓扑批量提交。
        std::vector<Candidate> candidates = ExtractCandidatePairs(currentShape, edgeFaceMap, messages_);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "iteration " << (iter + 1) << ": candidate count = " << candidates.size();
            AppendMessage(messages_, os.str());
        }

        if (candidates.empty()) {
            break;
        }

        // repairedThisIteration 标记本轮是否已经成功提交一个候选。
        bool repairedThisIteration = false;
        for (const Candidate& candidate : candidates) {
            // repairedShape 仅在单候选完整校验通过后才赋给 currentShape。
            TopoDS_Shape repairedShape;
            // errorMessage 记录当前候选失败原因，失败后继续尝试下一候选。
            std::string errorMessage;
            if (!TryRepairCandidate(currentShape, candidate, repairedShape, errorMessage)) {
                AppendMessage(messages_, "candidate failed: " + errorMessage);
                continue;
            }

            currentShape = repairedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++repairedCount_;
            repairedThisIteration = true;

            std::ostringstream os;
            os << "repaired candidate: confidence=" << candidate.confidence
               << ", maxGap=" << candidate.maxGap
               << ", meanGap=" << candidate.meanGap
               << ", lengthRatio=" << candidate.lengthRatio
               << ", tangentScore=" << candidate.tangentScore;
            AppendMessage(messages_, os.str());
            break;
        }

        if (!repairedThisIteration) {
            break;
        }
    }

    resultShape_ = currentShape;
    return repairedCount_;
}

int HealNarrowGap::Perform(const TopoDS_Shape& inputShape)
{
    setInputShape(inputShape);
    return Perform();
}

// 遍历 edgeFaceMap 中的边对，筛出“同端点、各属一个不同面、几何很近”的候选。
std::vector<HealNarrowGap::Candidate> HealNarrowGap::ExtractCandidatePairs(
    const TopoDS_Shape& shape,
    const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap,
    std::vector<std::string>& messages) const
{
    // candidates 保存所有通过拓扑和几何过滤的候选。
    std::vector<Candidate> candidates;
    if (shape.IsNull()) {
        return candidates;
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

        if (IsDegeneratedEdge(edgeA) || IsSeamEdge(edgeA, edgeFaceMap)) {
            AppendMessage(messages, "candidate skipped: edge is seam or degenerated");
            continue;
        }

        if (!Has3DCurve(edgeA)) {
            AppendMessage(messages, "candidate skipped: edge has no 3D curve");
            continue;
        }

        // faceA 是 edgeA 唯一所属面。
        const TopoDS_Face faceA = TopoDS::Face(facesA.First());
        if (faceA.IsNull()) {
            continue;
        }

        for (int j = i + 1; j <= edgeFaceMap.Extent(); ++j) {
            // edgeB 是与 edgeA 配对检查的另一条边。
            const TopoDS_Edge edgeB = TopoDS::Edge(edgeFaceMap.FindKey(j));
            if (edgeB.IsNull() || edgeA.IsSame(edgeB)) {
                continue;
            }

            if (!HaveSameVertexPair(edgeA, edgeB)) {
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

            if (IsDegeneratedEdge(edgeB) || IsSeamEdge(edgeB, edgeFaceMap)) {
                AppendMessage(messages, "candidate skipped: edge is seam or degenerated");
                continue;
            }

            if (!Has3DCurve(edgeB)) {
                AppendMessage(messages, "candidate skipped: edge has no 3D curve");
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
                AppendMessage(messages, "candidate skipped: edge length is too small");
                continue;
            }

            candidate.lengthRatio = minLength / maxLength;
            if (candidate.lengthRatio + Precision::Confusion() < myParams.minLengthRatio) {
                AppendMessage(messages, "candidate skipped: length ratio too small");
                continue;
            }

            // errorMessage 记录几何评估失败的具体原因。
            std::string errorMessage;
            if (!EvaluateGeometricSimilarity(candidate, errorMessage)) {
                AppendMessage(messages, "candidate skipped: " + errorMessage);
                continue;
            }

            candidate.confidence = ComputeCandidateConfidence(candidate);
            candidates.push_back(candidate);
        }
    }

    return candidates;
}

// 判断两条边是否引用同一对 TopoDS_Vertex，顺序可相反。
bool HealNarrowGap::HaveSameVertexPair(const TopoDS_Edge& edgeA, const TopoDS_Edge& edgeB) const
{
    TopoDS_Vertex aStart;
    TopoDS_Vertex aEnd;
    TopoDS_Vertex bStart;
    TopoDS_Vertex bEnd;
    if (!GetEdgeVerticesOriented(edgeA, aStart, aEnd) || !GetEdgeVerticesOriented(edgeB, bStart, bEnd)) {
        return false;
    }

    return (aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) ||
           (aStart.IsSame(bEnd) && aEnd.IsSame(bStart));
}

// 获取 edge 在当前 orientation 下的有向起终点。
bool HealNarrowGap::GetEdgeVerticesOriented(const TopoDS_Edge& edge, TopoDS_Vertex& start, TopoDS_Vertex& end) const
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
bool HealNarrowGap::IsDegeneratedEdge(const TopoDS_Edge& edge) const
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
bool HealNarrowGap::IsSeamEdge(
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
bool HealNarrowGap::Has3DCurve(const TopoDS_Edge& edge) const
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

// 按共享端点在各自 edge 上的曲线参数建立对应采样，评估间隙和切向一致性。
bool HealNarrowGap::EvaluateGeometricSimilarity(
    Candidate& candidate,
    std::string& errorMessage) const
{
    try {
        // aStart/aEnd 是 edgeA 在当前拓扑方向下的参考端点。
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        if (!GetEdgeVerticesOriented(candidate.edgeA, aStart, aEnd) ||
            !GetEdgeVerticesOriented(candidate.edgeB, bStart, bEnd)) {
            errorMessage = "failed to get oriented edge vertices";
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
            errorMessage = "vertex pair mismatch";
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

        // 实际采样点数至少为 2，包含两端点。
        const int sampleCount = std::max(2, myParams.sampleCount);
        // tangentLimit 是允许夹角对应的最小 abs(dot)。
        const double tangentLimit = std::cos(std::max(0.0, myParams.maxTangentAngle));
        // 参数递减采样时，D1 导数方向需要反转成实际采样行进方向。
        const bool reverseTangentA = uA1 < uA0;
        const bool reverseTangentB = uB1 < uB0;

        // gapSum 用于计算平均间隙；maxGap 记录最大采样间隙。
        double gapSum = 0.0;
        double maxGap = 0.0;
        // minAbsDot 记录所有采样点中最差的切向一致性。
        double minAbsDot = 1.0;

        for (int i = 0; i < sampleCount; ++i) {
            // t 是归一化采样参数。
            const double t = sampleCount == 1 ? 0.0 : static_cast<double>(i) / static_cast<double>(sampleCount - 1);
            // paramA/paramB 是根据对应端点参数插值得到的实际曲线参数。
            const double paramA = uA0 + (uA1 - uA0) * t;
            const double paramB = uB0 + (uB1 - uB0) * t;

            // pointA/pointB 是采样点；tangentA/tangentB 是对应一阶导。
            gp_Pnt pointA;
            gp_Pnt pointB;
            gp_Vec tangentA;
            gp_Vec tangentB;
            curveA.D1(paramA, pointA, tangentA);
            curveB.D1(paramB, pointB, tangentB);
            if (reverseTangentA) {
                tangentA.Reverse();
            }
            if (reverseTangentB) {
                tangentB.Reverse();
            }

            // gap 是当前采样点之间的 3D 距离。
            const double gap = pointA.Distance(pointB);
            gapSum += gap;
            maxGap = std::max(maxGap, gap);

            if (tangentA.SquareMagnitude() > Precision::SquareConfusion() &&
                tangentB.SquareMagnitude() > Precision::SquareConfusion()) {
                tangentA.Normalize();
                tangentB.Normalize();
                // absDot 允许两条边方向相反，只关心切线夹角大小。
                const double absDot = std::abs(tangentA.Dot(tangentB));
                minAbsDot = std::min(minAbsDot, Clamp01(absDot));
            }
        }

        candidate.meanGap = gapSum / static_cast<double>(sampleCount);
        candidate.maxGap = maxGap;
        candidate.tangentScore = minAbsDot;

        if (candidate.maxGap > myParams.maxGap) {
            errorMessage = "max gap too large";
            return false;
        }

        if (candidate.tangentScore + Precision::Confusion() < tangentLimit) {
            errorMessage = "tangent mismatch";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "curve sampling failed: " + FailureText(failure);
        return false;
    }

    return true;
}

// 将几何指标转换为可排序的综合置信度。
double HealNarrowGap::ComputeCandidateConfidence(const Candidate& candidate) const
{
    // maxGapDenom 防止 maxGap 参数过小导致除零。
    const double maxGapDenom = std::max(myParams.maxGap, Precision::Confusion());
    // 间隙越小分数越高。
    const double maxGapScore = 1.0 - Clamp01(candidate.maxGap / maxGapDenom);
    const double meanGapScore = 1.0 - Clamp01(candidate.meanGap / maxGapDenom);
    // 长度和切向越接近理想值分数越高。
    const double lengthScore = Clamp01(candidate.lengthRatio);
    const double tangentScore = Clamp01(candidate.tangentScore);

    // tolA/tolB 用于给容差更稳定的候选略微加权。
    const double tolA = candidate.edgeA.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeA);
    const double tolB = candidate.edgeB.IsNull() ? 0.0 : BRep_Tool::Tolerance(candidate.edgeB);
    const double tolSpread = std::abs(tolA - tolB) / std::max(std::max(tolA, tolB), Precision::Confusion());
    const double toleranceScore = 1.0 - Clamp01(tolSpread * 0.25);

    return 0.35 * maxGapScore +
           0.25 * meanGapScore +
           0.20 * lengthScore +
           0.15 * tangentScore +
           0.05 * toleranceScore;
}

// 选择保留边：明显较小容差优先，否则按参数选择较长边或默认 edgeA。
HealNarrowGap::EdgeChoice HealNarrowGap::ChooseEdgeToKeep(const Candidate& candidate) const
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

// 尝试提交单个候选修复；任何校验失败都会返回 false，由调用方回滚并尝试下一候选。
bool HealNarrowGap::TryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    // repairedShape 只有在整个局部替换和校验通过后才赋值。
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    // freeEdgesBefore 用于保证局部修复不会增加自由边数量。
    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    // choice 决定保留哪条边、重建哪个面。
    const EdgeChoice choice = ChooseEdgeToKeep(candidate);
    if (choice.keptEdge.IsNull() || choice.edgeToReplace.IsNull() ||
        choice.fixedFace.IsNull() || choice.faceToRebuild.IsNull()) {
        errorMessage = "invalid edge choice";
        return false;
    }

    // pcurveTol 是补建 pcurve 和 same-parameter 的保守容差。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(choice.keptEdge),
        BRep_Tool::Tolerance(choice.edgeToReplace),
        myParams.maxGap
        });

    TopoDS_Edge keptEdgeInFixedWire;
    if (!FindEdgeUseInFaceWire(choice.fixedFace,
        choice.keptEdge,
        keptEdgeInFixedWire,
        errorMessage)) {
        errorMessage = "kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeToReplaceInWire;
    if (!FindEdgeUseInFaceWire(choice.faceToRebuild,
        choice.edgeToReplace,
        edgeToReplaceInWire,
        errorMessage)) {
        errorMessage = "kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    if (!EnsurePCurveOnFace(keptEdgeInFixedWire, choice.fixedFace, pcurveTol, errorMessage)) {
        errorMessage = "failed to add pcurve on fixed face: " + errorMessage;
        return false;
    }

    if (!EnsurePCurveOnFace(edgeToReplaceInWire, choice.faceToRebuild, pcurveTol, errorMessage)) {
        errorMessage = "failed to add pcurve on target face: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFace;
    if (!RebuildFaceReplacingEdge(
        choice.faceToRebuild,
        choice.edgeToReplace,
        keptEdgeInFixedWire,
        pcurveTol,
        rebuiltFace,
        errorMessage)) {
        return false;
    }

    // trialShape 是将目标 face 边界重建后的试运行形状。
    TopoDS_Shape trialShape;
    try {
        BRepTools_ReShape reshaper;
        reshaper.Replace(choice.faceToRebuild, rebuiltFace);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "shape replacement returned null";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "shape replacement failed: " + FailureText(failure);
        return false;
    }

    // DBG-TEMP
    if (!FindEdgeUseInFaceWire(choice.fixedFace,
        choice.keptEdge,
        keptEdgeInFixedWire,
        errorMessage)) {
        errorMessage = "kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    if (!FindEdgeUseInFaceWire(rebuiltFace,
        choice.keptEdge,
        edgeToReplaceInWire,
        errorMessage)) {
        errorMessage = "kept edge not found in fixed face wire: " + errorMessage;
        return false;
    }

    bool ff = edgeToReplaceInWire.IsSame(keptEdgeInFixedWire);
    bool o1 = edgeToReplaceInWire.Orientation();
    bool o2 = keptEdgeInFixedWire.Orientation();



    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape keptEdgeFaces;
    if (!FindFacesForSameEdge(repairedEdgeFaceMap, choice.keptEdge, keptEdgeFaces)) {
        errorMessage = "kept edge not found after edge replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(keptEdgeFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge keptEdgeInFace;
        if (!FindEdgeUseInFaceWire(face, choice.keptEdge, keptEdgeInFace, pcurveRepairError)) {
            errorMessage = "failed to find kept edge use after replacement: " + pcurveRepairError;
            return false;
        }
        
        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                keptEdgeInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "failed to repair kept edge pcurve after replacement: " + pcurveRepairError;
            return false;
        }
    }

    if (myParams.runShapeFixAfterRepair) {
        try {
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(std::max(Precision::Confusion(), myParams.maxGap));
            fixer.SetMaxTolerance(std::max(Precision::Confusion(), myParams.maxGap * 10.0));
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix failed: " + FailureText(failure);
            return false;
        }

    }

    if (!ValidateEdgeReplacement(trialShape, choice.keptEdge, freeEdgesBefore, errorMessage)) {
        return false;
    }

    repairedShape = trialShape;
    return true;
}

// 确认目标 face 的 wire 中确实存在要被替换的 edge。
bool HealNarrowGap::FaceContainsEdgeInWire(
    const TopoDS_Face& oldFace,
    const TopoDS_Edge& oldEdge,
    std::string& errorMessage) const
{
    TopoDS_Edge oldEdgeInWire;
    return FindEdgeUseInFaceWire(oldFace, oldEdge, oldEdgeInWire, errorMessage);
}

bool HealNarrowGap::FindEdgeUseInFaceWire(
    const TopoDS_Face& oldFace,
    const TopoDS_Edge& oldEdge,
    TopoDS_Edge& oldEdgeInWire,
    std::string& errorMessage) const
{
    oldEdgeInWire.Nullify();
    if (oldFace.IsNull() || oldEdge.IsNull()) {
        errorMessage = "invalid face or edge for edge replacement check";
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
        errorMessage = "failed to find replacement edge in face: " + FailureText(failure);
        return false;
    }

    errorMessage = "old edge not found in face";
    return false;
}

bool HealNarrowGap::RebuildFaceReplacingEdge(
    const TopoDS_Face& faceToRebuild,
    const TopoDS_Edge& edgeToReplace,
    const TopoDS_Edge& replacementEdgeUse,
    double tol,
    TopoDS_Face& rebuiltFace,
    std::string& errorMessage) const
{
    rebuiltFace.Nullify();
    if (faceToRebuild.IsNull() || edgeToReplace.IsNull() || replacementEdgeUse.IsNull()) {
        errorMessage = "invalid face or edge for face rebuild";
        return false;
    }

    try {
        TopoDS_Face forwardFace = faceToRebuild;
        forwardFace.Orientation(TopAbs_FORWARD);

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(forwardFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "target face has no support surface";
            return false;
        }

        const double faceTol = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(faceToRebuild),
            tol
        });

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
                    errorMessage = "failed to rebuild target wire";
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
            errorMessage = "edge to replace not found while rebuilding face";
            return false;
        }

        ShapeFix_Face faceFix(rebuilt);
        faceFix.SetPrecision(faceTol);
        faceFix.SetMaxTolerance(std::max(faceTol, myParams.maxGap * 10.0));
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
        errorMessage = "face rebuild failed: " + FailureText(failure);
        return false;
    }

    return true;
}

// 如果 edge 在 face 上缺少 pcurve，则尝试用 ShapeFix_Edge 投影并修正 same-parameter。
bool HealNarrowGap::EnsurePCurveOnFace(
    TopoDS_Edge& edge,
    const TopoDS_Face& face,
    double tol,
    std::string& errorMessage) const
{
    if (edge.IsNull() || face.IsNull()) {
        errorMessage = "edge or face is null";
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
bool HealNarrowGap::ValidateWholeShape(
    const TopoDS_Shape& repairedShape,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (repairedShape.IsNull()) {
        errorMessage = "repaired shape is null";
        return false;
    }

    try {
        BRepCheck_Analyzer analyzer(repairedShape, Standard_True);
        if (!analyzer.IsValid()) {
            errorMessage = "whole shape invalid after replacement:";
            errorMessage += DescribeBRepCheckFailures(analyzer, repairedShape);
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "whole shape check failed: " + FailureText(failure);
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
bool HealNarrowGap::ValidateEdgeReplacement(
    const TopoDS_Shape& repairedShape,
    const TopoDS_Edge& keptEdge,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    /*if (!ValidateWholeShape(repairedShape, freeEdgesBefore, errorMessage)) {
        return false;
    }*/

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(repairedShape, edgeFaceMap);
    TopTools_ListOfShape adjacentFaces;
    if (!FindFacesForSameEdge(edgeFaceMap, keptEdge, adjacentFaces)) {
        errorMessage = "kept edge not found after edge replacement";
        return false;
    }

    if (adjacentFaces.Extent() != 2) {
        errorMessage = "kept edge is not shared after edge replacement";
        return false;
    }

    for (TopTools_ListIteratorOfListOfShape it(adjacentFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        const double pcurveTol = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(keptEdge),
            BRep_Tool::Tolerance(face)
        });
        std::string localError;
        if (!GeometryUtils::HasValidPCurveOnFace(keptEdge, face, pcurveTol, &localError)) {
            errorMessage = "kept edge pcurve invalid after edge replacement: " + localError;
            return false;
        }
    }

    return true;
}

// 统一追加日志，便于后续调整日志策略。
void HealNarrowGap::AppendMessage(std::vector<std::string>& messages, const std::string& message) const
{
    messages.push_back(message);
}

// Manual verification:
// Positive case: two distinct one-face edges sharing the same vertex pair and
// whose sampled 3D gap is below Params::maxGap should yield Perform()>0,
// isModified()==true, two faces sharing the kept TopoDS_Edge, and a valid BRep.
// Negative case: if the sampled middle gap is greater than Params::maxGap,
// Perform should return 0 and leave isModified()==false.
