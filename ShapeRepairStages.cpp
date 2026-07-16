#include "ShapeRepairStages.h"

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
