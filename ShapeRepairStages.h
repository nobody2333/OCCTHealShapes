#pragma once

#include "HealClosedNarrowGap.h"
#include "HealNarrowGap.h"
#include "HealSimpleHole.h"
#include "ShapeRepairPipeline.h"
#include "SewShells.h"
#include "MergeDuplicateVertex.h"

// 将当前已有修复类适配为通用流水线阶段。
// 新模块可以在这里增加便捷工厂，也可以在调用处直接构造自定义 Stage。
namespace ShapeRepairStages
{
    ShapeRepairPipeline::Stage MakeMergeDuplicateVertexStage(
        const MergeDuplicateVertex::Params& params = MergeDuplicateVertex::Params(),
        bool enabled = true);

    ShapeRepairPipeline::Stage MakeConservativeSewShellsStage(
        double tolerance = 1.0e-4, bool enabled = true);

    ShapeRepairPipeline::Stage MakeHealNarrowGapStage(
        const HealNarrowGap::Params& params = HealNarrowGap::Params(),
        bool enabled = true);

    ShapeRepairPipeline::Stage MakeHealClosedNarrowGapStage(
        const HealClosedNarrowGap::Params& params = HealClosedNarrowGap::Params(),
        bool enabled = true);

    ShapeRepairPipeline::Stage MakeHealSimpleHoleStage(
        const HealSimpleHole::Params& params = HealSimpleHole::Params(),
        bool enabled = true);
}
