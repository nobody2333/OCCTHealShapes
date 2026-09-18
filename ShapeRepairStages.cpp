#include "ShapeRepairStages.h"

ShapeRepairPipeline::Stage ShapeRepairStages::MakeShapePreProcessStage(
    const ShapePreProcess::Params& params, bool enabled)
{
    ShapeRepairPipeline::Stage stage;
    stage.name = "ShapePreProcess";
    stage.enabled = enabled;
    stage.execute = [params](const TopoDS_Shape& inputShape) {
        ShapePreProcess preProcess(params);
        ShapeRepairPipeline::StageExecutionResult result;
        result.repairCount = preProcess.perform(inputShape);
        result.resultShape = preProcess.getResultShape();
        result.modified = preProcess.isModified();
        result.messages = preProcess.getMessages();
        return result;
    };
    return stage;
}

ShapeRepairPipeline::Stage ShapeRepairStages::MakeConservativeSewShellsStage(
    double tolerance, bool enabled)
{
    ShapeRepairPipeline::Stage stage;
    stage.name = "ConservativeSewShells";
    stage.enabled = enabled;
    stage.execute = [tolerance](const TopoDS_Shape& inputShape) {
        SewShells sewer(SewShells::ConservativeParams(tolerance));
        ShapeRepairPipeline::StageExecutionResult result;
        result.repairCount = sewer.perform(inputShape);
        result.resultShape = sewer.getResultShape();
        result.modified = sewer.isModified();
        result.messages = sewer.getMessages();
        return result;
    };
    return stage;
}

ShapeRepairPipeline::Stage ShapeRepairStages::MakeHealNarrowGapStage(
    const HealNarrowGap::Params& params,
    bool enabled)
{
    ShapeRepairPipeline::Stage stage;
    stage.name = "HealNarrowGap";
    stage.enabled = enabled;
    stage.execute = [params](const TopoDS_Shape& inputShape) {
        HealNarrowGap healer(params);
        ShapeRepairPipeline::StageExecutionResult result;
        result.repairCount = healer.perform(inputShape);
        result.resultShape = healer.getResultShape();
        result.modified = healer.isModified();
        result.messages = healer.getMessages();
        return result;
    };
    return stage;
}

ShapeRepairPipeline::Stage ShapeRepairStages::MakeHealClosedNarrowGapStage(
    const HealClosedNarrowGap::Params& params,
    bool enabled)
{
    ShapeRepairPipeline::Stage stage;
    stage.name = "HealClosedNarrowGap";
    stage.enabled = enabled;
    stage.execute = [params](const TopoDS_Shape& inputShape) {
        HealClosedNarrowGap healer(params);
        ShapeRepairPipeline::StageExecutionResult result;
        result.repairCount = healer.perform(inputShape);
        result.resultShape = healer.getResultShape();
        result.modified = healer.isModified();
        result.messages = healer.getMessages();
        return result;
    };
    return stage;
}

ShapeRepairPipeline::Stage ShapeRepairStages::MakeHealSimpleHoleStage(
    const HealSimpleHole::Params& params,
    bool enabled)
{
    ShapeRepairPipeline::Stage stage;
    stage.name = "HealSimpleHole";
    stage.enabled = enabled;
    stage.execute = [params](const TopoDS_Shape& inputShape) {
        HealSimpleHole healer(params);
        ShapeRepairPipeline::StageExecutionResult result;
        result.repairCount = healer.perform(inputShape);
        result.resultShape = healer.getResultShape();
        result.modified = healer.isModified();
        result.messages = healer.getMessages();
        return result;
    };
    return stage;
}
