#include "SewShells.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
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

namespace {
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

SewShells::SewShells(const Params& params)
    : params_(params)
{
    params_.sewingTolerance = std::max(params_.sewingTolerance, Precision::Confusion());
    params_.minTolerance = std::max(params_.minTolerance, Precision::Confusion());
    params_.maxTolerance = std::max(params_.maxTolerance, params_.sewingTolerance);
    params_.minShellCountPerGroup = std::max(2, params_.minShellCountPerGroup);
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

    for (TopoDS_Iterator childIterator(sourceShape); childIterator.More(); childIterator.Next()) {
        const TopoDS_Shape childShape = childIterator.Value();
        if (childShape.IsNull()) {
            continue;
        }

        if (childShape.ShapeType() == TopAbs_SHELL) {
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
        for (const TopoDS_Shell& shell : shellGroup) {
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
            shellFixer.Perform();

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
        for (TopoDS_Iterator childIterator(shape); childIterator.More(); childIterator.Next()) {
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
