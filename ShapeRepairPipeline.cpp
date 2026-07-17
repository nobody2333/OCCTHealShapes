#include "ShapeRepairPipeline.h"

#include <algorithm>
#include <exception>
#include <sstream>
#include <utility>

#include <Standard_Failure.hxx>

#include "RepairDiagnostics.h"

namespace
{
    std::string standardFailureMessage(const Standard_Failure& failure)
    {
        const char* message = failure.GetMessageString();
        return message != nullptr ? message : "unknown OCCT failure";
    }

    bool collectFailures(
        const TopoDS_Shape& shape,
        RepairDiagnostics::FailureList& failures,
        std::string& errorMessage)
    {
        errorMessage.clear();
        failures = RepairDiagnostics::CollectBRepCheckFailures(shape, &errorMessage);
        return errorMessage.empty();
    }

    std::string describeFailureCountIncrease(
        const RepairDiagnostics::FailureList& existingFailures,
        const RepairDiagnostics::FailureList& resultFailures)
    {
        const RepairDiagnostics::FailureList newFailures =
            RepairDiagnostics::FindBRepCheckFailureCountIncrease(
                existingFailures,
                resultFailures);
        if (newFailures.empty()) {
            return {};
        }

        return "new BRepCheck failures detected: "
            + RepairDiagnostics::DescribeBRepCheckFailures(newFailures);
    }
}

ShapeRepairPipeline::ShapeRepairPipeline(const Params& params)
    : params_(params)
{
}

void ShapeRepairPipeline::addStage(const Stage& stage)
{
    stages_.push_back(stage);
}

void ShapeRepairPipeline::addStage(Stage&& stage)
{
    stages_.push_back(std::move(stage));
}

void ShapeRepairPipeline::clearStages()
{
    stages_.clear();
}

ShapeRepairPipeline::Result ShapeRepairPipeline::perform(const TopoDS_Shape& inputShape) const
{
    Result pipelineResult;
    pipelineResult.resultShape = inputShape;

    if (inputShape.IsNull()) {
        pipelineResult.messages.push_back("Pipeline input shape is null");
        return pipelineResult;
    }

    RepairDiagnostics::FailureList originalFailures;
    std::string originalFailureError;
    const bool hasOriginalFailureBaseline = collectFailures(
        inputShape,
        originalFailures,
        originalFailureError);
    if (!hasOriginalFailureBaseline) {
        pipelineResult.messages.push_back(
            "Cannot establish pipeline input BRepCheck baseline: "
            + originalFailureError);
        if (params_.validateFinalResult || params_.validateAfterEachStage) {
            return pipelineResult;
        }
    }

    TopoDS_Shape workingShape = inputShape;
    for (const Stage& stage : stages_) {
        StageReport report;
        report.name = stage.name.empty() ? "UnnamedStage" : stage.name;
        report.enabled = stage.enabled;
        report.inputShapeType = workingShape.ShapeType();
        report.outputShapeType = workingShape.ShapeType();

        if (!stage.enabled) {
            report.messages.push_back("Stage disabled by pipeline configuration");
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        if (!stage.execute) {
            report.rejectionReason = "stage executor is empty";
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        RepairDiagnostics::FailureList stageInputFailures;
        std::string stageInputFailureError;
        const bool hasStageFailureBaseline = !params_.validateAfterEachStage
            || collectFailures(workingShape, stageInputFailures, stageInputFailureError);
        if (params_.validateAfterEachStage && !hasStageFailureBaseline) {
            report.rejectionReason =
                "failed to collect stage input BRepCheck failures: "
                + stageInputFailureError;
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        StageExecutionResult executionResult;
        report.executed = true;
        try {
            executionResult = stage.execute(workingShape);
        }
        catch (const Standard_Failure& failure) {
            report.rejectionReason =
                "OCCT failure while executing stage: " + standardFailureMessage(failure);
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }
        catch (const std::exception& exception) {
            report.rejectionReason =
                "standard exception while executing stage: " + std::string(exception.what());
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }
        catch (...) {
            report.rejectionReason = "unknown exception while executing stage";
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        report.modified = executionResult.modified;
        report.repairCount = std::max(0, executionResult.repairCount);
        report.messages = std::move(executionResult.messages);
        if (!executionResult.resultShape.IsNull()) {
            report.outputShapeType = executionResult.resultShape.ShapeType();
        }

        // 未修改阶段是正常结果，后续阶段继续使用当前 workingShape。
        if (!executionResult.modified) {
            report.accepted = true;
            if (executionResult.repairCount != 0) {
                report.messages.push_back(
                    "Warning: stage reported repairs but modified=false; repair count was not committed");
            }
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        if (executionResult.resultShape.IsNull()) {
            report.rejectionReason = "modified stage returned a null result shape";
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        if (params_.rejectShapeTypeChange
            && executionResult.resultShape.ShapeType() != workingShape.ShapeType()) {
            std::ostringstream reason;
            reason << "stage changed top-level shape type from "
                   << static_cast<int>(workingShape.ShapeType())
                   << " to "
                   << static_cast<int>(executionResult.resultShape.ShapeType());
            report.rejectionReason = reason.str();
            pipelineResult.stageReports.push_back(std::move(report));
            continue;
        }

        if (params_.validateAfterEachStage) {
            RepairDiagnostics::FailureList stageResultFailures;
            std::string stageResultFailureError;
            if (!collectFailures(
                    executionResult.resultShape,
                    stageResultFailures,
                    stageResultFailureError)) {
                report.rejectionReason =
                    "failed to collect stage result BRepCheck failures: "
                    + stageResultFailureError;
                pipelineResult.stageReports.push_back(std::move(report));
                continue;
            }

            if (hasStageFailureBaseline) {
                report.rejectionReason = describeFailureCountIncrease(
                    stageInputFailures,
                    stageResultFailures);
                if (!report.rejectionReason.empty()) {
                    pipelineResult.stageReports.push_back(std::move(report));
                    continue;
                }
            }
        }

        // 只有通过流水线公共校验的阶段结果才会成为下一阶段输入。
        workingShape = executionResult.resultShape;
        report.accepted = true;
        pipelineResult.modified = true;
        pipelineResult.totalRepairCount += report.repairCount;
        pipelineResult.stageReports.push_back(std::move(report));
    }

    pipelineResult.resultShape = workingShape;

    if (params_.validateFinalResult
        && pipelineResult.modified
        && hasOriginalFailureBaseline) {
        RepairDiagnostics::FailureList finalFailures;
        std::string finalFailureError;
        if (!collectFailures(workingShape, finalFailures, finalFailureError)) {
            pipelineResult.messages.push_back(
                "Final pipeline validation failed: " + finalFailureError);
            pipelineResult.resultShape = inputShape;
            pipelineResult.modified = false;
            pipelineResult.rolledBack = true;
            pipelineResult.totalRepairCount = 0;
        }
        else {
            const std::string newFailureDescription = describeFailureCountIncrease(
                originalFailures,
                finalFailures);
            if (!newFailureDescription.empty()) {
                pipelineResult.messages.push_back(
                    "Final pipeline result rejected: " + newFailureDescription);
                pipelineResult.resultShape = inputShape;
                pipelineResult.modified = false;
                pipelineResult.rolledBack = true;
                pipelineResult.totalRepairCount = 0;
            }
        }
    }

    return pipelineResult;
}
