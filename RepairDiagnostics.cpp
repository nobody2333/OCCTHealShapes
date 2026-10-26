#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <map>
#include <sstream>

#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

namespace {
bool containsStatus(
    const std::vector<BRepCheck_Status>& statuses,
    BRepCheck_Status status)
{
    return std::find(statuses.begin(), statuses.end(), status) != statuses.end();
}

bool statusesAreSubset(
    const std::vector<BRepCheck_Status>& statuses,
    const std::vector<BRepCheck_Status>& existingStatuses)
{
    for (const BRepCheck_Status status : statuses) {
        if (!containsStatus(existingStatuses, status)) {
            return false;
        }
    }
    return true;
}

std::string replaceAll(
    std::string text,
    const std::string& from,
    const std::string& to)
{
    if (from.empty()) {
        return text;
    }

    std::string::size_type pos = 0;
    while ((pos = text.find(from, pos)) != std::string::npos) {
        text.replace(pos, from.length(), to);
        pos += to.length();
    }
    return text;
}
} // namespace

const char* RepairDiagnostics::BRepCheckStatusName(BRepCheck_Status status)
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

RepairDiagnostics::FailureList RepairDiagnostics::CollectBRepCheckFailures(
    const TopoDS_Shape& shape,
    std::string* errorMessage)
{
    FailureList failures;
    if (errorMessage != nullptr) {
        errorMessage->clear();
    }

    if (shape.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "shape is null";
        }
        return failures;
    }

    try {
        // 使用 OCCT 自带的并行检查；几何检查和精度模式保持原值，
        // 错误仍在检查完成后按原有拓扑遍历顺序串行收集。
        BRepCheck_Analyzer analyzer(shape, Standard_True, Standard_True);
        TopTools_IndexedMapOfShape allShapes;
        TopExp::MapShapes(shape, allShapes);

        for (int i = 1; i <= allShapes.Extent(); ++i) {
            const TopoDS_Shape& subShape = allShapes.FindKey(i);
            Handle(BRepCheck_Result) result = analyzer.Result(subShape);
            if (result.IsNull()) {
                continue;
            }

            BRepCheckFailureInfo failure;
            failure.shape = subShape;
            failure.type = subShape.ShapeType();

            const BRepCheck_ListOfStatus& statuses = result->Status();
            for (BRepCheck_ListIteratorOfListOfStatus it(statuses); it.More(); it.Next()) {
                const BRepCheck_Status status = it.Value();
                if (status == BRepCheck_NoError || containsStatus(failure.statuses, status)) {
                    continue;
                }
                failure.statuses.push_back(status);
            }

            if (!failure.statuses.empty()) {
                failures.push_back(failure);
            }
        }
    }
    catch (Standard_Failure& failure) {
        if (errorMessage != nullptr) {
            *errorMessage = "BRepCheck failure collection failed: " + FailureText(failure);
        }
    }

    return failures;
}

bool RepairDiagnostics::IsPreExistingBRepCheckFailure(
    const BRepCheckFailureInfo& failure,
    const FailureList& existingFailures)
{
    for (const BRepCheckFailureInfo& existingFailure : existingFailures) {
        if (failure.type != existingFailure.type || !failure.shape.IsSame(existingFailure.shape)) {
            continue;
        }

        if (statusesAreSubset(failure.statuses, existingFailure.statuses)) {
            return true;
        }
    }
    return false;
}

RepairDiagnostics::FailureList RepairDiagnostics::FindNewBRepCheckFailures(
    const FailureList& existingFailures,
    const FailureList& resultFailures)
{
    FailureList newFailures;
    for (const BRepCheckFailureInfo& failure : resultFailures) {
        if (!IsPreExistingBRepCheckFailure(failure, existingFailures)) {
            newFailures.push_back(failure);
        }
    }
    return newFailures;
}

RepairDiagnostics::FailureList RepairDiagnostics::FindBRepCheckFailureCountIncrease(
    const FailureList& existingFailures,
    const FailureList& resultFailures)
{
    using FailureKey = std::pair<int, int>;
    std::map<FailureKey, int> remainingExistingCounts;
    for (const BRepCheckFailureInfo& failure : existingFailures) {
        for (const BRepCheck_Status status : failure.statuses) {
            ++remainingExistingCounts[FailureKey(
                static_cast<int>(failure.type),
                static_cast<int>(status))];
        }
    }

    FailureList increasedFailures;
    for (const BRepCheckFailureInfo& failure : resultFailures) {
        BRepCheckFailureInfo increasedFailure;
        increasedFailure.shape = failure.shape;
        increasedFailure.type = failure.type;

        for (const BRepCheck_Status status : failure.statuses) {
            const FailureKey key(
                static_cast<int>(failure.type),
                static_cast<int>(status));
            auto existingIt = remainingExistingCounts.find(key);
            if (existingIt != remainingExistingCounts.end() && existingIt->second > 0) {
                --existingIt->second;
            }
            else {
                increasedFailure.statuses.push_back(status);
            }
        }

        if (!increasedFailure.statuses.empty()) {
            increasedFailures.push_back(increasedFailure);
        }
    }

    return increasedFailures;
}

std::string RepairDiagnostics::DescribeBRepCheckFailures(
    const FailureList& failures,
    int maxReported)
{
    std::ostringstream os;
    int reported = 0;
    maxReported = std::max(0, maxReported);

    for (const BRepCheckFailureInfo& failure : failures) {
        if (reported >= maxReported) {
            os << " ...";
            break;
        }

        os << " [type=" << static_cast<int>(failure.type) << " status=";
        for (std::size_t i = 0; i < failure.statuses.size(); ++i) {
            if (i > 0) {
                os << ",";
            }
            os << BRepCheckStatusName(failure.statuses[i]);
        }
        os << "]";
        ++reported;
    }
    return os.str();
}

std::string RepairDiagnostics::DescribeBRepCheckFailures(
    BRepCheck_Analyzer& analyzer,
    const TopoDS_Shape& shape,
    int maxReported)
{
    std::ostringstream os;
    TopTools_IndexedMapOfShape allShapes;
    TopExp::MapShapes(shape, allShapes);

    int reported = 0;
    maxReported = std::max(0, maxReported);
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

    if (reported == maxReported && allShapes.Extent() > reported) {
        os << " ...";
    }
    return os.str();
}

std::string RepairDiagnostics::MakeFailureReasonReadable(const std::string& reason)
{
    if (reason.empty()) {
        return "no detail";
    }

    std::string readable = replaceAll(reason, "; ", "\n");
    readable = replaceAll(readable, ": ", ":\n  ");
    return readable;
}

std::string RepairDiagnostics::IndentMultiline(
    const std::string& text,
    const std::string& indent)
{
    std::ostringstream os;
    bool atLineStart = true;
    for (char ch : text) {
        if (atLineStart) {
            os << indent;
        }
        os << ch;
        atLineStart = ch == '\n';
    }
    return os.str();
}

std::string RepairDiagnostics::StrategyFailureBlock(
    const std::string& strategyName,
    const std::string& reason)
{
    std::ostringstream os;
    os << "[" << strategyName << "]\n"
       << IndentMultiline(MakeFailureReasonReadable(reason), "  ");
    return os.str();
}

void RepairDiagnostics::AppendIterationHeader(
    std::vector<std::string>& messages,
    int iteration,
    bool addLeadingBlankLine)
{
    if (addLeadingBlankLine) {
        messages.emplace_back();
    }

    messages.push_back("=== Iteration " + std::to_string(iteration) + " ===");
}

void RepairDiagnostics::AppendShapeCheckSummary(
    std::vector<std::string>& messages,
    const std::string& shapeLabel,
    const TopoDS_Shape& shape)
{
    std::string errorMessage;
    const FailureList failures = CollectBRepCheckFailures(shape, &errorMessage);
    if (!errorMessage.empty()) {
        messages.push_back(shapeLabel + " BRepCheck failed: " + errorMessage);
        return;
    }

    if (failures.empty()) {
        messages.push_back(shapeLabel + " has no BRepCheck failures");
        return;
    }

    messages.push_back(
        shapeLabel + " has existing BRepCheck failures:"
        + DescribeBRepCheckFailures(failures));
}

std::string RepairDiagnostics::FailureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
}

std::string RepairDiagnostics::DescribeEdgeFaceCounts(
    const TopoDS_Shape& shape,
    const std::vector<TopoDS_Edge>& edges)
{
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);

    std::ostringstream os;
    os << "edgeCount=" << edges.size() << ", adjacentFaceCounts=[";
    for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
        if (edgeIndex > 0) {
            os << ", ";
        }

        int faceCount = 0;
        if (!edges[edgeIndex].IsNull() && edgeFaceMap.Contains(edges[edgeIndex])) {
            faceCount = edgeFaceMap.FindFromKey(edges[edgeIndex]).Extent();
        }
        os << edgeIndex << ":" << faceCount;
    }
    os << "]";
    return os.str();
}
