#include "SewShells.h"
#include "TopologyUtils.h"
#include "RepairDiagnostics.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopoDS_Vertex.hxx>

namespace {
double subshapeTolerance(const TopoDS_Shape& shape)
{
    switch (shape.ShapeType()) {
    case TopAbs_VERTEX: return BRep_Tool::Tolerance(TopoDS::Vertex(shape));
    case TopAbs_EDGE: return BRep_Tool::Tolerance(TopoDS::Edge(shape));
    case TopAbs_FACE: return BRep_Tool::Tolerance(TopoDS::Face(shape));
    default: return 0.0;
    }
}

std::string failureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
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

struct BRepCheckFailureInfo {
    TopoDS_Shape shape;
    TopAbs_ShapeEnum type = TopAbs_SHAPE;
    std::vector<BRepCheck_Status> statuses;
};

bool containsStatus(const std::vector<BRepCheck_Status>& statuses, BRepCheck_Status status)
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

std::vector<BRepCheckFailureInfo> collectBRepCheckFailures(
    const TopoDS_Shape& shape,
    std::string* errorMessage = nullptr)
{
    std::vector<BRepCheckFailureInfo> failures;
    if (shape.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "shape is null";
        }
        return failures;
    }

    try {
        BRepCheck_Analyzer analyzer(shape, Standard_True);
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
            *errorMessage = "BRepCheck failure collection failed: " + failureText(failure);
        }
    }

    return failures;
}

bool isPreExistingBRepCheckFailure(
    const BRepCheckFailureInfo& failure,
    const std::vector<BRepCheckFailureInfo>& existingFailures)
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

std::string describeBRepCheckFailureList(const std::vector<BRepCheckFailureInfo>& failures)
{
    std::ostringstream os;
    int reported = 0;
    constexpr int maxReported = 30;
    for (const BRepCheckFailureInfo& failure : failures) {
        if (reported >= maxReported) {
            os << " ...";
            break;
        }

        os << " [type=" << TopologyUtils::ShapeTypeName(failure.type) << " status=";
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
} // namespace

SewShells::Params SewShells::ConservativeParams(double tolerance)
{
    Params params;
    params.conservative = true;
    params.sewingTolerance = tolerance;
    return params;
}

SewShells::SewShells(const Params& params)
    : params_(params)
{
    if (!std::isfinite(params_.sewingTolerance)) params_.sewingTolerance = 1.0e-4;
    if (!std::isfinite(params_.minTolerance)) params_.minTolerance = Precision::Confusion();
    if (!std::isfinite(params_.maxTolerance)) params_.maxTolerance = params_.sewingTolerance;
    params_.sewingTolerance = std::max(params_.sewingTolerance, Precision::Confusion());
    params_.minTolerance = std::max(params_.minTolerance, Precision::Confusion());
    params_.maxTolerance = std::max(params_.maxTolerance, params_.sewingTolerance);
    params_.minShellCountPerGroup = std::max(2, params_.minShellCountPerGroup);
    params_.minTolerance = std::min(params_.minTolerance, params_.sewingTolerance);
    if (params_.conservative) {
        params_.maxTolerance = params_.sewingTolerance;
        params_.minTolerance = Precision::Confusion();
        params_.localTolerancesMode = false;
        params_.nonManifoldMode = false;
        params_.floatingEdgesMode = false;
        params_.faceMode = true;
        params_.sameParameterMode = true;
        params_.runShapeFixAfterSewing = false;
        params_.discardLooseTopLevelSubshapesAfterSewing = false;
        params_.rejectFreeEdgeIncrease = true;
        params_.rejectInvalidResult = true;
    }
}

void SewShells::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = inputShape_;
    modified_ = false;
    sewedGroupCount_ = 0;
    messages_.clear();
}

int SewShells::perform()
{
    resultShape_ = inputShape_;
    modified_ = false;
    sewedGroupCount_ = 0;
    messages_.clear();

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    {
        std::ostringstream os;
        os << "Start SewShells: inputType=" << TopologyUtils::ShapeTypeName(inputShape_)
           << ", conservative=" << (params_.conservative ? "true" : "false")
           << ", sewingTolerance=" << params_.sewingTolerance
           << ", minTolerance=" << params_.minTolerance
           << ", maxTolerance=" << params_.maxTolerance
           << ", minShellCountPerGroup=" << params_.minShellCountPerGroup
           << ", faceMode=" << (params_.faceMode ? "true" : "false")
           << ", sameParameterMode=" << (params_.sameParameterMode ? "true" : "false")
           << ", floatingEdgesMode=" << (params_.floatingEdgesMode ? "true" : "false")
           << ", localTolerancesMode=" << (params_.localTolerancesMode ? "true" : "false")
           << ", nonManifoldMode=" << (params_.nonManifoldMode ? "true" : "false")
           << ", fixShellOrientation=" << (params_.fixShellOrientation ? "true" : "false")
           << ", requireSingleShell=" << (params_.requireSingleShell ? "true" : "false")
           << ", runShapeFixAfterSewing=" << (params_.runShapeFixAfterSewing ? "true" : "false")
           << ", discardLooseTopLevelSubshapesAfterSewing="
           << (params_.discardLooseTopLevelSubshapesAfterSewing ? "true" : "false");
        appendMessage(os.str());
    }

    TopoDS_Shape sewedShape;
    int sewedCount = 0;
    if (!sewShapeRecursively(inputShape_, sewedShape, sewedCount) || sewedShape.IsNull()) {
        appendMessage("SewShells stopped: recursive sewing failed");
        resultShape_ = inputShape_;
        return 0;
    }

    sewedGroupCount_ = sewedCount;
    modified_ = sewedGroupCount_ > 0;
    if (modified_ && params_.runShapeFixAfterSewing) {
        std::string shapeFixError;
        if (!applyShapeFixAfterSewing(sewedShape, "whole sewed shape", shapeFixError)) {
            appendMessage("Warning: " + shapeFixError);
        }
    }

    resultShape_ = sewedShape;

    {
        std::ostringstream os;
        os << "Finish SewShells: sewedGroupCount=" << sewedGroupCount_
           << ", modified=" << (modified_ ? "true" : "false")
           << ", resultType=" << TopologyUtils::ShapeTypeName(resultShape_);
        appendMessage(os.str());
    }

    return sewedGroupCount_;
}

int SewShells::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

bool SewShells::sewShapeRecursively(
    const TopoDS_Shape& sourceShape,
    TopoDS_Shape& resultShape,
    int& sewedGroupCount)
{
    resultShape = sourceShape;
    if (sourceShape.IsNull()) {
        return true;
    }

    // 当前实现只在 compound 层级内寻找相同父label shell。
    // solid 内部 shell 不会被递归修改，避免破坏已有实体边界。
    if (sourceShape.ShapeType() != TopAbs_COMPOUND) {
        return true;
    }

    std::vector<TopoDS_Shape> rebuiltChildren;
    std::vector<TopoDS_Shell> directShells;
    bool childWasModified = false;

    // 子项保持局部位置/方向，重建后再恢复父级，避免重复施加变换。
    for (TopoDS_Iterator childIterator(sourceShape, Standard_False, Standard_False); childIterator.More(); childIterator.Next()) {
        const TopoDS_Shape childShape = childIterator.Value();
        if (childShape.IsNull()) {
            continue;
        }

        if (childShape.ShapeType() == TopAbs_SHELL) {
            if (params_.conservative && BRep_Tool::IsClosed(childShape)) {
                rebuiltChildren.push_back(childShape);
                appendMessage("Conservative sewing: closed shell preserved");
                continue;
            }
            directShells.push_back(TopoDS::Shell(childShape));
            continue;
        }

        TopoDS_Shape rebuiltChild = childShape;
        const int beforeChildSewedCount = sewedGroupCount;
        if (!sewShapeRecursively(childShape, rebuiltChild, sewedGroupCount)) {
            rebuiltChild = childShape;
        }

        if (sewedGroupCount != beforeChildSewedCount) {
            childWasModified = true;
        }
        rebuiltChildren.push_back(rebuiltChild);
    }

    bool directShellsWereSewed = false;
    TopoDS_Shape directSewedShape;
    if (static_cast<int>(directShells.size()) >= params_.minShellCountPerGroup) {
        std::string errorMessage;
        if (sewShellGroup(directShells, directSewedShape, errorMessage)) {
            directShellsWereSewed = true;
            ++sewedGroupCount;
        }
        else {
            appendMessage("Shell group sewing failed: " + errorMessage);
        }
    }

    if (!directShellsWereSewed && !childWasModified) {
        if (params_.conservative && static_cast<int>(directShells.size()) < params_.minShellCountPerGroup) {
            appendMessage("Conservative sewing skipped: fewer than two direct open shells in compound");
        }
        resultShape = sourceShape;
        return true;
    }

    BRep_Builder builder;
    TopoDS_Compound rebuiltCompound;
    builder.MakeCompound(rebuiltCompound);

    for (const TopoDS_Shape& childShape : rebuiltChildren) {
        if (!childShape.IsNull()) {
            builder.Add(rebuiltCompound, childShape);
        }
    }

    if (directShellsWereSewed) {
        if (!directSewedShape.IsNull()) {
            if (directSewedShape.ShapeType() == TopAbs_COMPOUND) {
                bool addedChild = false;
                for (TopoDS_Iterator sewedChildIterator(directSewedShape);
                     sewedChildIterator.More();
                     sewedChildIterator.Next()) {
                    const TopoDS_Shape sewedChild = sewedChildIterator.Value();
                    if (!sewedChild.IsNull()) {
                        builder.Add(rebuiltCompound, sewedChild);
                        addedChild = true;
                    }
                }
                if (!addedChild) {
                    builder.Add(rebuiltCompound, directSewedShape);
                }
            }
            else {
                builder.Add(rebuiltCompound, directSewedShape);
            }
        }
    }
    else {
        for (const TopoDS_Shell& shell : directShells) {
            if (!shell.IsNull()) {
                builder.Add(rebuiltCompound, shell);
            }
        }
    }

    rebuiltCompound.Orientation(sourceShape.Orientation());
    rebuiltCompound.Location(sourceShape.Location());
    resultShape = rebuiltCompound;
    return true;
}

bool SewShells::sewShellGroup(
    const std::vector<TopoDS_Shell>& shellGroup,
    TopoDS_Shape& sewedShape,
    std::string& errorMessage)
{
    if (params_.conservative) {
        return sewCrossShellBoundaries(shellGroup, sewedShape, errorMessage);
    }
    sewedShape.Nullify();
    if (static_cast<int>(shellGroup.size()) < params_.minShellCountPerGroup) {
        errorMessage = "not enough shells in group";
        return false;
    }

    const TopoDS_Shape beforeShape = makeShellGroupCompound(shellGroup);
    if (beforeShape.IsNull()) {
        errorMessage = "shell group compound is null";
        return false;
    }

    try {
        // 一次复制整个组，保留组内共享拓扑，失败不会污染输入的曲线和容差。
        TopoDS_Shape trialShape = beforeShape;
        if (params_.conservative) {
            BRepBuilderAPI_Copy trialCopy(beforeShape, Standard_True, Standard_False);
            trialShape = trialCopy.Shape();
            if (trialShape.IsNull()) {
                errorMessage = "shell group trial copy is null";
                return false;
            }
        }
        // 在 sewing 修改副本之前保存逐子形状容差。只有身份未变且容差未增的
        // 原有子形状可以保留超预算容差，新建/替换子形状仍严格受预算限制。
        TopTools_IndexedMapOfShape trialSubshapes;
        std::vector<double> trialTolerances;
        if (params_.conservative) {
            TopExp::MapShapes(trialShape, trialSubshapes);
            for (int i = 1; i <= trialSubshapes.Extent(); ++i) {
                trialTolerances.push_back(subshapeTolerance(trialSubshapes(i)));
            }
        }
        BRepBuilderAPI_Sewing sewing(
            params_.sewingTolerance,
            params_.faceMode,
            Standard_True,
            params_.floatingEdgesMode,
            params_.nonManifoldMode);

        sewing.SetTolerance(params_.sewingTolerance);
        sewing.SetMinTolerance(params_.minTolerance);
        sewing.SetMaxTolerance(params_.maxTolerance);
        sewing.SetFaceMode(params_.faceMode);
        sewing.SetSameParameterMode(params_.sameParameterMode);
        sewing.SetFloatingEdgesMode(params_.floatingEdgesMode);
        sewing.SetLocalTolerancesMode(params_.localTolerancesMode);
        sewing.SetNonManifoldMode(params_.nonManifoldMode);

        int inputFaceCount = 0;
        for (TopExp_Explorer shellIt(trialShape, TopAbs_SHELL); shellIt.More(); shellIt.Next()) {
            const TopoDS_Shell shell = TopoDS::Shell(shellIt.Current());
            if (shell.IsNull()) {
                continue;
            }
            sewing.Add(shell);

            for (TopExp_Explorer faceExplorer(shell, TopAbs_FACE);
                faceExplorer.More();
                faceExplorer.Next()) {
                ++inputFaceCount;
            }
        }

        if (inputFaceCount == 0) {
            errorMessage = "Shell group has no faces";
            return false;
        }

        sewing.Perform();
        TopoDS_Shape rawSewedShape = sewing.SewedShape();
        if (params_.conservative && (sewing.NbDeletedFaces() != 0 || sewing.NbMultipleEdges() != 0)) {
            errorMessage = "conservative sewing rejected deleted faces or non-manifold edges";
            return false;
        }
        if (rawSewedShape.IsNull()) {
            errorMessage = "BRepBuilderAPI_Sewing produced a null shape";
            return false;
        }

        {
            std::ostringstream os;
            os << "Raw sewing result: resultType=" << TopologyUtils::ShapeTypeName(rawSewedShape)
               << ", shellCount=" << countShells(rawSewedShape)
               << ", " << describeDirectChildren(rawSewedShape)
               << ", freeEdges=" << sewing.NbFreeEdges()
               << ", multipleEdges=" << sewing.NbMultipleEdges()
               << ", contiguousEdges=" << sewing.NbContigousEdges()
               << ", deletedFaces=" << sewing.NbDeletedFaces();
            appendMessage(os.str());
        }

        std::string orientationWarning;
        if (params_.fixShellOrientation) {
            rawSewedShape = fixShellOrientationsInShape(rawSewedShape, orientationWarning);
            if (!orientationWarning.empty()) {
                appendMessage("Warning: " + orientationWarning);
            }

            std::ostringstream os;
            os << "ShapeFix_Shell result: resultType=" << TopologyUtils::ShapeTypeName(rawSewedShape)
               << ", shellCount=" << countShells(rawSewedShape)
               << ", " << describeDirectChildren(rawSewedShape);
            appendMessage(os.str());
        }

        if (!handleLooseTopLevelSubshapes(rawSewedShape, "raw sewing result", errorMessage)) {
            return false;
        }

        TopoDS_Shape acceptedShape;
        int outputShellCount = 0;
        if (!acceptSewedResult(
                rawSewedShape,
                static_cast<int>(shellGroup.size()),
                acceptedShape,
                outputShellCount,
                errorMessage)) {
            return false;
        }

        if (params_.runShapeFixAfterSewing) {
            std::string shapeFixError;
            if (!applyShapeFixAfterSewing(acceptedShape, "sewed shell group", shapeFixError)) {
                errorMessage = shapeFixError;
                return false;
            }

            if (!handleLooseTopLevelSubshapes(
                    acceptedShape,
                    "ShapeFix_Shape result after sewing",
                    errorMessage)) {
                return false;
            }

            TopoDS_Shape shapeFixAcceptedShape;
            if (!acceptSewedResult(
                    acceptedShape,
                    static_cast<int>(shellGroup.size()),
                    shapeFixAcceptedShape,
                    outputShellCount,
                    errorMessage)) {
                errorMessage = "ShapeFix_Shape result is not acceptable after sewing: " + errorMessage;
                return false;
            }
            acceptedShape = shapeFixAcceptedShape;
        }

        if (params_.conservative) {
            TopTools_IndexedMapOfShape outputSubshapes;
            TopExp::MapShapes(acceptedShape, outputSubshapes);
            int preservedOverBudget = 0;
            for (int i = 1; i <= outputSubshapes.Extent(); ++i) {
                const TopoDS_Shape& subshape = outputSubshapes(i);
                const double tolerance = subshapeTolerance(subshape);
                const int sourceIndex = trialSubshapes.FindIndex(subshape);
                const double originalTolerance = sourceIndex > 0
                    ? trialTolerances[static_cast<std::size_t>(sourceIndex - 1)] : 0.0;
                const bool unchangedTolerance = sourceIndex > 0
                    && std::isfinite(originalTolerance) && tolerance <= originalTolerance;
                if (!std::isfinite(tolerance)
                    || (tolerance > params_.maxTolerance && !unchangedTolerance)) {
                    std::ostringstream reason;
                    reason << "conservative sewing new/increased tolerance exceeds budget: type="
                        << TopologyUtils::ShapeTypeName(subshape) << ", tolerance=" << tolerance
                        << ", sameTrialSubshape=" << (sourceIndex > 0)
                        << ", originalTolerance=" << originalTolerance
                        << ", budget=" << params_.maxTolerance;
                    errorMessage = reason.str();
                    return false;
                }
                if (tolerance > params_.maxTolerance) ++preservedOverBudget;
            }
            appendMessage("Tolerance validation: preserved unchanged over-budget subshapes="
                + std::to_string(preservedOverBudget));
        }
        std::string validationWarning;
        if (!validateSewingResult(beforeShape, acceptedShape, errorMessage, validationWarning)) {
            return false;
        }
        if (!validationWarning.empty()) {
            appendMessage("Warning: " + validationWarning);
        }

        {
            std::ostringstream os;
            os << "Sewed shell group: inputShells=" << shellGroup.size()
               << ", inputFaces=" << inputFaceCount
               << ", rawOutputShells=" << countShells(rawSewedShape)
               << ", acceptedOutputShells=" << outputShellCount
               << ", acceptedDirectChildren=" << countDirectChildren(acceptedShape)
               << ", requireSingleShell=" << (params_.requireSingleShell ? "true" : "false")
               << ", resultType=" << TopologyUtils::ShapeTypeName(acceptedShape)
               << ", " << describeDirectChildren(acceptedShape)
               << ", freeEdges=" << sewing.NbFreeEdges()
               << ", multipleEdges=" << sewing.NbMultipleEdges()
               << ", contiguousEdges=" << sewing.NbContigousEdges()
               << ", deletedFaces=" << sewing.NbDeletedFaces();
            appendMessage(os.str());
        }

        sewedShape = acceptedShape;
        return true;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "OCCT failure while sewing shells: " + failureText(failure);
        return false;
    }
}

bool SewShells::sewCrossShellBoundaries(
    const std::vector<TopoDS_Shell>& shellGroup,
    TopoDS_Shape& sewedShape,
    std::string& errorMessage)
{
    sewedShape.Nullify();
    const TopoDS_Shape original = makeShellGroupCompound(shellGroup);
    TopoDS_Shape current = original;
    // 每个下标永久代表输入组中的同一子形状；提交后只更新其像，不重新判定归属。
    TopTools_IndexedMapOfShape sources;
    TopExp::MapShapes(original, sources);
    std::vector<TopoDS_Shape> images;
    for (int i = 1; i <= sources.Extent(); ++i) images.push_back(sources(i));
    std::vector<std::vector<int>> shellMembers(shellGroup.size());
    for (std::size_t s = 0; s < shellGroup.size(); ++s) {
        TopTools_IndexedMapOfShape members;
        TopExp::MapShapes(shellGroup[s], members);
        for (int i = 1; i <= members.Extent(); ++i) {
            const auto type = members(i).ShapeType();
            const int id = sources.FindIndex(members(i)) - 1;
            if (type == TopAbs_EDGE || type == TopAbs_VERTEX || type == TopAbs_FACE)
                shellMembers[s].push_back(id);
        }
    }
    struct Boundary { int edge; int face; std::size_t shell; };
    std::vector<Boundary> boundaries;
    TopTools_IndexedDataMapOfShapeListOfShape originalAncestors;
    TopExp::MapShapesAndAncestors(original, TopAbs_EDGE, TopAbs_FACE, originalAncestors);
    for (std::size_t s = 0; s < shellGroup.size(); ++s) {
        for (const int id : shellMembers[s]) {
            if (images[id].ShapeType() != TopAbs_EDGE) continue;
            const auto edge = TopoDS::Edge(images[id]);
            const auto& faces = originalAncestors.FindFromKey(edge);
            Standard_Real first = 0.0, last = 0.0;
            if (faces.Extent() != 1 || BRep_Tool::Degenerated(edge)
                || BRep_Tool::Curve(edge, first, last).IsNull()
                || !std::isfinite(first) || !std::isfinite(last) || last <= first) continue;
            TopoDS_Vertex a, b;
            TopExp::Vertices(edge, a, b);
            if (a.IsNull() || b.IsNull() || a.IsSame(b)
                || BRep_Tool::IsClosed(edge, TopoDS::Face(faces.First()))) continue;
            boundaries.push_back({id, sources.FindIndex(faces.First()) - 1, s});
        }
    }
    struct Pair { std::size_t a; std::size_t b; double distance; };
    std::vector<Pair> pairs;
    for (std::size_t a = 0; a < boundaries.size(); ++a) {
        TopoDS_Vertex a0, a1;
        TopExp::Vertices(TopoDS::Edge(images[boundaries[a].edge]), a0, a1);
        for (std::size_t b = a + 1; b < boundaries.size(); ++b) {
            if (boundaries[a].shell == boundaries[b].shell
                || boundaries[a].edge == boundaries[b].edge) continue;
            TopoDS_Vertex b0, b1;
            TopExp::Vertices(TopoDS::Edge(images[boundaries[b].edge]), b0, b1);
            const double distance = std::min(
                std::max(BRep_Tool::Pnt(a0).Distance(BRep_Tool::Pnt(b0)),
                         BRep_Tool::Pnt(a1).Distance(BRep_Tool::Pnt(b1))),
                std::max(BRep_Tool::Pnt(a0).Distance(BRep_Tool::Pnt(b1)),
                         BRep_Tool::Pnt(a1).Distance(BRep_Tool::Pnt(b0))));
            if (distance <= params_.sewingTolerance) pairs.push_back({a, b, distance});
        }
    }
    std::stable_sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
        return a.distance < b.distance;
    });
    appendMessage("Cross-shell local sewing: originalShells=" + std::to_string(shellGroup.size())
        + ", boundaryEdges=" + std::to_string(boundaries.size())
        + ", endpointCandidates=" + std::to_string(pairs.size()));
    int committed = 0;
    for (const Pair& pair : pairs) {
        const Boundary& a = boundaries[pair.a];
        const Boundary& b = boundaries[pair.b];
        if (images[a.edge].IsSame(images[b.edge])) continue;
        std::string reason;
        try {
            TopTools_IndexedDataMapOfShapeListOfShape ancestors;
            TopExp::MapShapesAndAncestors(current, TopAbs_EDGE, TopAbs_FACE, ancestors);
            if (!ancestors.Contains(images[a.edge]) || !ancestors.Contains(images[b.edge])
                || ancestors.FindFromKey(images[a.edge]).Extent() != 1
                || ancestors.FindFromKey(images[b.edge]).Extent() != 1) continue;
            BRepBuilderAPI_Copy copy(current, Standard_True, Standard_False);
            std::vector<TopoDS_Shape> trialImages(images.size());
            for (std::size_t i = 0; i < images.size(); ++i) {
                const auto type = images[i].ShapeType();
                if (type == TopAbs_EDGE || type == TopAbs_VERTEX || type == TopAbs_FACE)
                    trialImages[i] = copy.ModifiedShape(images[i]);
            }
            TopTools_IndexedMapOfShape trialShapes;
            TopExp::MapShapes(copy.Shape(), trialShapes);
            std::vector<double> tolerances;
            for (int i = 1; i <= trialShapes.Extent(); ++i)
                tolerances.push_back(subshapeTolerance(trialShapes(i)));

            // 只 Add 两个候选面，不 Load 整组；关闭 cutting，禁止一对多切边。
            BRepBuilderAPI_Sewing sewing(params_.sewingTolerance, Standard_True,
                Standard_True, Standard_False, Standard_False);
            sewing.SetMinTolerance(params_.minTolerance);
            sewing.SetMaxTolerance(params_.maxTolerance);
            sewing.SetSameParameterMode(params_.sameParameterMode);
            sewing.SetLocalTolerancesMode(Standard_False);
            sewing.Add(trialImages[a.face]);
            sewing.Add(trialImages[b.face]);
            sewing.Perform();
            if (sewing.NbDeletedFaces() || sewing.NbMultipleEdges())
                reason = "deleted face or non-manifold edge";
            const auto context = sewing.GetContext();
            for (auto& shape : trialImages) {
                if (shape.IsNull()) continue;
                const auto type = shape.ShapeType();
                shape = context->Apply(shape);
                if (shape.IsNull() || shape.ShapeType() != type) {
                    reason = "deleted/split subshape";
                    break;
                }
            }
            if (reason.empty() && !trialImages[a.edge].IsSame(trialImages[b.edge]))
                reason = "target edges were not joined";
            // 核心约束：每个原始 shell 内，原本不同的边/点/面必须仍然不同。
            for (std::size_t s = 0; reason.empty() && s < shellMembers.size(); ++s) {
                TopTools_IndexedMapOfShape unique;
                for (const int id : shellMembers[s]) {
                    if (unique.Contains(trialImages[id])) {
                        reason = "would merge distinct subshapes within original shell "
                            + std::to_string(s + 1);
                        break;
                    }
                    unique.Add(trialImages[id]);
                }
            }
            // 局部面内的其它边也可能被 OCCT 匹配；只允许原本的自由边参与连接。
            TopTools_IndexedMapOfShape mappedEdges;
            std::vector<int> firstSources;
            for (std::size_t i = 0; reason.empty() && i < images.size(); ++i) {
                if (images[i].ShapeType() != TopAbs_EDGE) continue;
                const int index = mappedEdges.FindIndex(trialImages[i]);
                if (index == 0) {
                    mappedEdges.Add(trialImages[i]);
                    firstSources.push_back(static_cast<int>(i));
                } else {
                    const auto& previous = images[firstSources[index - 1]];
                    if (!previous.IsSame(images[i])
                        && (ancestors.FindFromKey(previous).Extent() != 1
                            || ancestors.FindFromKey(images[i]).Extent() != 1))
                        reason = "would merge an already connected edge";
                }
            }
            if (reason.empty()) {
                // 仅按共享边合并壳容器，不进行任何基于距离的全局 sewing。
                std::vector<std::size_t> components(shellGroup.size());
                for (std::size_t s = 0; s < components.size(); ++s) components[s] = s;
                TopTools_IndexedMapOfShape seenEdges;
                std::vector<std::size_t> edgeOwner;
                for (std::size_t s = 0; s < shellMembers.size(); ++s) {
                    for (const int id : shellMembers[s]) {
                        if (trialImages[id].ShapeType() != TopAbs_EDGE) continue;
                        const int index = seenEdges.FindIndex(trialImages[id]);
                        if (!index) {
                            seenEdges.Add(trialImages[id]);
                            edgeOwner.push_back(s);
                        } else {
                            const auto from = components[s];
                            const auto to = components[edgeOwner[index - 1]];
                            for (auto& component : components) if (component == from) component = to;
                        }
                    }
                }
                BRep_Builder builder;
                TopoDS_Compound rebuilt;
                builder.MakeCompound(rebuilt);
                for (std::size_t s = 0; s < components.size(); ++s) {
                    if (components[s] != s) continue;
                    TopoDS_Shell shell;
                    builder.MakeShell(shell);
                    for (std::size_t t = 0; t < components.size(); ++t) {
                        if (components[t] != s) continue;
                        // 从当前源面恢复累计方向；历史表的共享子形状方向不是 face use。
                        for (TopExp_Explorer f(shellGroup[t], TopAbs_FACE); f.More(); f.Next()) {
                            const int id = sources.FindIndex(f.Current()) - 1;
                            builder.Add(shell, trialImages[id].Oriented(f.Current().Orientation()));
                        }
                    }
                    builder.Add(rebuilt, shell);
                }
                TopoDS_Shape result = rebuilt;
                std::string warning;
                if (params_.fixShellOrientation) result = fixShellOrientationsInShape(result, warning);
                TopTools_IndexedMapOfShape output;
                TopExp::MapShapes(result, output);
                for (int i = 1; reason.empty() && i <= output.Extent(); ++i) {
                    const double tolerance = subshapeTolerance(output(i));
                    const int source = trialShapes.FindIndex(output(i));
                    if (!std::isfinite(tolerance) || (tolerance > params_.maxTolerance
                        && (source == 0 || tolerance > tolerances[source - 1])))
                        reason = "new/increased tolerance exceeds budget";
                }
                // 方向整理只能改变 use 方向，不能使历史表指向已丢失的子形状。
                for (const auto& shape : trialImages)
                    if (reason.empty() && !shape.IsNull() && !output.Contains(shape))
                        reason = "orientation repair changed tracked topology";
                if (reason.empty() && validateSewingResult(current, result, reason, warning)) {
                    current = result;
                    for (std::size_t i = 0; i < images.size(); ++i)
                        if (!trialImages[i].IsNull()) images[i] = trialImages[i];
                    ++committed;
                }
            }
        } catch (const Standard_Failure& failure) {
            reason = failureText(failure);
        }
        appendMessage("Local sewing shells " + std::to_string(a.shell + 1) + "/"
            + std::to_string(b.shell + 1) + ", edges " + std::to_string(a.edge + 1)
            + "/" + std::to_string(b.edge + 1) + ": "
            + (reason.empty() ? "committed" : "rejected: " + reason));
    }
    if (!committed) {
        errorMessage = "no admissible cross-shell local sewing";
        return false;
    }
    int outputShells = 0;
    TopoDS_Shape accepted;
    std::string warning;
    if (!acceptSewedResult(current, static_cast<int>(shellGroup.size()), accepted,
            outputShells, errorMessage)
        || !validateSewingResult(original, accepted, errorMessage, warning)) return false;
    sewedShape = accepted;
    appendMessage("Cross-shell local sewing committed pairs=" + std::to_string(committed)
        + ", outputShells=" + std::to_string(outputShells));
    return true;
}

TopoDS_Shape SewShells::makeShellGroupCompound(const std::vector<TopoDS_Shell>& shellGroup) const
{
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    for (const TopoDS_Shell& shell : shellGroup) {
        if (!shell.IsNull()) {
            builder.Add(compound, shell);
        }
    }

    return compound;
}

TopoDS_Shape SewShells::fixShellOrientationsInShape(
    const TopoDS_Shape& shape,
    std::string& warningMessage) const
{
    warningMessage.clear();
    if (shape.IsNull()) {
        return shape;
    }

    try {
        if (shape.ShapeType() == TopAbs_SHELL) {
            ShapeFix_Shell shellFixer(TopoDS::Shell(shape));
            if (params_.conservative) {
                shellFixer.FixFaceOrientation(TopoDS::Shell(shape), Standard_True, Standard_False);
            }
            else {
                shellFixer.Perform();
            }

            const TopoDS_Shape fixedShape = shellFixer.Shape();
            return fixedShape.IsNull() ? shape : fixedShape;
        }

        if (shape.ShapeType() != TopAbs_COMPOUND) {
            return shape;
        }

        BRep_Builder builder;
        TopoDS_Compound fixedCompound;
        builder.MakeCompound(fixedCompound);

        bool childWasFixed = false;
        for (TopoDS_Iterator childIterator(shape, Standard_False, Standard_False); childIterator.More(); childIterator.Next()) {
            const TopoDS_Shape childShape = childIterator.Value();
            std::string childWarning;
            TopoDS_Shape fixedChild = fixShellOrientationsInShape(childShape, childWarning);
            if (!childWarning.empty()) {
                if (!warningMessage.empty()) {
                    warningMessage += "; ";
                }
                warningMessage += childWarning;
            }

            if (!fixedChild.IsNull()) {
                builder.Add(fixedCompound, fixedChild);
                if (!fixedChild.IsSame(childShape)) {
                    childWasFixed = true;
                }
            }
        }

        if (!childWasFixed) {
            return shape;
        }

        fixedCompound.Orientation(shape.Orientation());
        fixedCompound.Location(shape.Location());
        return fixedCompound;
    }
    catch (Standard_Failure& failure) {
        warningMessage = "ShapeFix_Shell failed after sewing: " + failureText(failure);
        return shape;
    }
}

bool SewShells::acceptSewedResult(
    const TopoDS_Shape& rawSewedShape,
    int inputShellCount,
    TopoDS_Shape& acceptedShape,
    int& outputShellCount,
    std::string& errorMessage) const
{
    acceptedShape.Nullify();
    outputShellCount = countShells(rawSewedShape);
    if (rawSewedShape.IsNull()) {
        errorMessage = "sewing result is null";
        return false;
    }

    if (params_.requireSingleShell) {
        TopoDS_Shell mergedShell;
        if (!extractSingleShell(rawSewedShape, mergedShell, errorMessage)) {
            return false;
        }

        acceptedShape = mergedShell;
        outputShellCount = 1;
        return true;
    }

    if (outputShellCount <= 0) {
        std::ostringstream os;
        os << "sewing result contains no shell: resultType="
           << TopologyUtils::ShapeTypeName(rawSewedShape);
        errorMessage = os.str();
        return false;
    }

    if (outputShellCount >= inputShellCount) {
        std::ostringstream os;
        os << "sewing did not reduce shell count: inputShells=" << inputShellCount
           << ", outputShells=" << outputShellCount
           << ", resultType=" << TopologyUtils::ShapeTypeName(rawSewedShape)
           << ", " << describeDirectChildren(rawSewedShape);
        errorMessage = os.str();
        return false;
    }


    acceptedShape = rawSewedShape;
    return true;
}

bool SewShells::applyShapeFixAfterSewing(
    TopoDS_Shape& shape,
    const std::string& context,
    std::string& errorMessage)
{
    if (shape.IsNull()) {
        errorMessage = "ShapeFix_Shape skipped for " + context + ": shape is null";
        return false;
    }

    try {
        ShapeFix_Shape fixer(shape);
        const double precision = std::max(Precision::Confusion(), params_.sewingTolerance);
        fixer.SetPrecision(precision);
        fixer.SetMaxTolerance(std::max({ precision, params_.maxTolerance, params_.sewingTolerance }));
        fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
        fixer.FixFaceTool()->FixMissingSeamMode() = 1;
        fixer.Perform();

        const TopoDS_Shape fixedShape = fixer.Shape();
        if (fixedShape.IsNull()) {
            errorMessage = "ShapeFix_Shape produced a null shape for " + context;
            return false;
        }

        const int shellCountBefore = countShells(shape);
        const int shellCountAfter = countShells(fixedShape);
        shape = fixedShape;

        std::ostringstream os;
        os << "ShapeFix_Shape applied to " << context
           << ": shellCountBefore=" << shellCountBefore
           << ", shellCountAfter=" << shellCountAfter
           << ", resultType=" << TopologyUtils::ShapeTypeName(shape);
        appendMessage(os.str());
        return true;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "ShapeFix_Shape failed for " + context + ": " + failureText(failure);
        return false;
    }
}

bool SewShells::handleLooseTopLevelSubshapes(
    TopoDS_Shape& shape,
    const std::string& context,
    std::string& errorMessage)
{
    if (shape.IsNull() || !hasLooseTopLevelSubshape(shape)) {
        return true;
    }

    const std::string beforeDescription = describeDirectChildren(shape);
    if (!params_.discardLooseTopLevelSubshapesAfterSewing) {
        errorMessage = "sewing result contains loose top-level subshape in "
            + context + ": " + beforeDescription;
        return false;
    }

    if (shape.ShapeType() != TopAbs_COMPOUND) {
        errorMessage = "cannot discard loose top-level subshape because the whole "
            + context + " is " + TopologyUtils::ShapeTypeName(shape);
        return false;
    }

    BRep_Builder builder;
    TopoDS_Compound cleanedCompound;
    builder.MakeCompound(cleanedCompound);

    int keptCount = 0;
    int discardedCount = 0;
    for (TopoDS_Iterator childIterator(shape); childIterator.More(); childIterator.Next()) {
        const TopoDS_Shape childShape = childIterator.Value();
        if (childShape.IsNull()) {
            continue;
        }

        if (isLooseTopLevelSubshape(childShape)) {
            ++discardedCount;
            continue;
        }

        builder.Add(cleanedCompound, childShape);
        ++keptCount;
    }

    if (keptCount == 0) {
        errorMessage = "discarding loose top-level subshapes would make "
            + context + " empty: " + beforeDescription;
        return false;
    }

    cleanedCompound.Orientation(shape.Orientation());
    cleanedCompound.Location(shape.Location());
    shape = cleanedCompound;

    std::ostringstream os;
    os << "Warning: discarded loose top-level subshapes after sewing"
       << ": context=" << context
       << ", discarded=" << discardedCount
       << ", before=" << beforeDescription
       << ", after=" << describeDirectChildren(shape);
    appendMessage(os.str());
    return true;
}

bool SewShells::isLooseTopLevelSubshape(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return false;
    }

    const TopAbs_ShapeEnum shapeType = shape.ShapeType();
    return shapeType == TopAbs_FACE
        || shapeType == TopAbs_WIRE
        || shapeType == TopAbs_EDGE
        || shapeType == TopAbs_VERTEX;
}

bool SewShells::hasLooseTopLevelSubshape(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return false;
    }

    if (shape.ShapeType() != TopAbs_COMPOUND) {
        return isLooseTopLevelSubshape(shape);
    }

    for (TopoDS_Iterator childIterator(shape); childIterator.More(); childIterator.Next()) {
        if (isLooseTopLevelSubshape(childIterator.Value())) {
            return true;
        }
    }
    return false;
}

bool SewShells::extractSingleShell(
    const TopoDS_Shape& shape,
    TopoDS_Shell& shell,
    std::string& errorMessage) const
{
    shell.Nullify();
    if (shape.IsNull()) {
        errorMessage = "sewing result is null";
        return false;
    }

    if (shape.ShapeType() == TopAbs_SHELL) {
        shell = TopoDS::Shell(shape);
        return !shell.IsNull();
    }

    TopTools_IndexedMapOfShape shellMap;
    TopExp::MapShapes(shape, TopAbs_SHELL, shellMap);
    if (shellMap.Extent() != 1) {
        std::ostringstream os;
        os << "sewing result did not merge into a single shell: resultType="
           << TopologyUtils::ShapeTypeName(shape)
           << ", shellCount=" << shellMap.Extent();
        errorMessage = os.str();
        return false;
    }

    shell = TopoDS::Shell(shellMap.FindKey(1));
    if (shell.IsNull()) {
        errorMessage = "single shell extracted from sewing result is null";
        return false;
    }

    return true;
}

bool SewShells::validateSewingResult(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& afterShape,
    std::string& errorMessage,
    std::string& warningMessage) const
{
    warningMessage.clear();
    const auto appendWarning = [&warningMessage](const std::string& message) {
        if (!warningMessage.empty()) {
            warningMessage += "; ";
        }
        warningMessage += message;
    };

    if (beforeShape.IsNull() || afterShape.IsNull()) {
        errorMessage = "before or after sewing shape is null";
        return false;
    }

    if (params_.conservative) {
        TopTools_IndexedMapOfShape beforeFaces, afterFaces;
        TopExp::MapShapes(beforeShape, TopAbs_FACE, beforeFaces);
        TopExp::MapShapes(afterShape, TopAbs_FACE, afterFaces);
        if (beforeFaces.Extent() != afterFaces.Extent()) {
            errorMessage = "conservative sewing changed face count";
            return false;
        }
        std::string diagnosticError;
        const auto before = RepairDiagnostics::CollectBRepCheckFailures(beforeShape, &diagnosticError);
        if (!diagnosticError.empty()) { errorMessage = diagnosticError; return false; }
        const auto after = RepairDiagnostics::CollectBRepCheckFailures(afterShape, &diagnosticError);
        if (!diagnosticError.empty()) { errorMessage = diagnosticError; return false; }
        const auto added = RepairDiagnostics::FindBRepCheckFailureCountIncrease(before, after);
        if (!added.empty()) {
            errorMessage = "conservative sewing added BRepCheck errors: "
                + RepairDiagnostics::DescribeBRepCheckFailures(added);
            return false;
        }
        const int beforeFree = TopologyUtils::CountFreeEdges(beforeShape);
        const int afterFree = TopologyUtils::CountFreeEdges(afterShape);
        if (afterFree >= beforeFree) {
            errorMessage = "conservative sewing did not reduce free edges";
            return false;
        }
        return true;
    }

    std::string existingError;
    const std::vector<BRepCheckFailureInfo> existingFailures =
        collectBRepCheckFailures(beforeShape, &existingError);
    if (!existingError.empty()) {
        errorMessage = "existing BRepCheck failure collection failed before sewing: " + existingError;
        return false;
    }

    if (!existingFailures.empty()) {
        appendWarning("existing BRepCheck failures before sewing:"
            + describeBRepCheckFailureList(existingFailures));
    }

    std::string afterError;
    const std::vector<BRepCheckFailureInfo> afterFailures =
        collectBRepCheckFailures(afterShape, &afterError);
    if (!afterError.empty()) {
        errorMessage = "sewed shape BRepCheck failed: " + afterError;
        return false;
    }

    std::vector<BRepCheckFailureInfo> newFailures;
    for (const BRepCheckFailureInfo& failure : afterFailures) {
        if (!isPreExistingBRepCheckFailure(failure, existingFailures)) {
            newFailures.push_back(failure);
        }
    }

    if (!newFailures.empty()) {
        std::string message = "new BRepCheck failures after sewing:";
        message += describeBRepCheckFailureList(newFailures);
        if (!existingFailures.empty()) {
            message += " existing failures:";
            message += describeBRepCheckFailureList(existingFailures);
        }

        if (params_.rejectInvalidResult) {
            errorMessage = message;
            return false;
        }
        appendWarning(message);
    }
    else if (!afterFailures.empty()) {
        appendWarning("sewed shape keeps only existing BRepCheck failures:"
            + describeBRepCheckFailureList(afterFailures));
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(beforeShape);
    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (freeEdgesAfter > freeEdgesBefore) {
        std::ostringstream os;
        os << "free edge count increased after sewing: before=" << freeEdgesBefore
           << ", after=" << freeEdgesAfter;
        if (params_.rejectFreeEdgeIncrease) {
            errorMessage = os.str();
            return false;
        }
        appendWarning(os.str());
    }

    return true;
}

int SewShells::countShells(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return 0;
    }

    if (shape.ShapeType() == TopAbs_SHELL) {
        return 1;
    }

    int count = 0;
    for (TopExp_Explorer shellExplorer(shape, TopAbs_SHELL);
         shellExplorer.More();
         shellExplorer.Next()) {
        ++count;
    }
    return count;
}

int SewShells::countDirectChildren(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return 0;
    }

    if (shape.ShapeType() != TopAbs_COMPOUND) {
        return 1;
    }

    int count = 0;
    for (TopoDS_Iterator childIterator(shape); childIterator.More(); childIterator.Next()) {
        if (!childIterator.Value().IsNull()) {
            ++count;
        }
    }
    return count;
}

std::string SewShells::describeDirectChildren(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return "directChildren=0 []";
    }

    if (shape.ShapeType() != TopAbs_COMPOUND) {
        std::ostringstream os;
        os << "directChildren=1 [" << TopologyUtils::ShapeTypeName(shape) << "=1]";
        return os.str();
    }

    int compoundCount = 0;
    int compSolidCount = 0;
    int solidCount = 0;
    int shellCount = 0;
    int faceCount = 0;
    int wireCount = 0;
    int edgeCount = 0;
    int vertexCount = 0;
    int unknownCount = 0;

    for (TopoDS_Iterator childIterator(shape); childIterator.More(); childIterator.Next()) {
        const TopoDS_Shape childShape = childIterator.Value();
        if (childShape.IsNull()) {
            continue;
        }

        switch (childShape.ShapeType()) {
        case TopAbs_COMPOUND: ++compoundCount; break;
        case TopAbs_COMPSOLID: ++compSolidCount; break;
        case TopAbs_SOLID: ++solidCount; break;
        case TopAbs_SHELL: ++shellCount; break;
        case TopAbs_FACE: ++faceCount; break;
        case TopAbs_WIRE: ++wireCount; break;
        case TopAbs_EDGE: ++edgeCount; break;
        case TopAbs_VERTEX: ++vertexCount; break;
        default: ++unknownCount; break;
        }
    }

    const int totalCount =
        compoundCount + compSolidCount + solidCount + shellCount + faceCount
        + wireCount + edgeCount + vertexCount + unknownCount;

    std::ostringstream os;
    os << "directChildren=" << totalCount << " [";
    bool first = true;
    const auto appendCount = [&os, &first](const char* name, int count) {
        if (count <= 0) {
            return;
        }
        if (!first) {
            os << ", ";
        }
        os << name << "=" << count;
        first = false;
    };

    appendCount("COMPOUND", compoundCount);
    appendCount("COMPSOLID", compSolidCount);
    appendCount("SOLID", solidCount);
    appendCount("SHELL", shellCount);
    appendCount("FACE", faceCount);
    appendCount("WIRE", wireCount);
    appendCount("EDGE", edgeCount);
    appendCount("VERTEX", vertexCount);
    appendCount("UNKNOWN", unknownCount);

    os << "]";
    return os.str();
}

void SewShells::appendMessage(const std::string& message)
{
    messages_.push_back(message);
}
