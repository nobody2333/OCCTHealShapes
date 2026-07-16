#pragma once

#include <functional>
#include <string>
#include <vector>

#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>

// ShapeRepairPipeline 负责按添加顺序串行执行多个 shape 修复阶段。
// 每个阶段只接收上一阶段已提交的结果，不负责 XCAF label 或 STEP 文件读写。
class ShapeRepairPipeline
{
public:
    struct Params
    {
        // 为 true 时，拒绝改变输入顶层 shape 类型的阶段结果。
        bool rejectShapeTypeChange = true;

        // 为 true 时，每个阶段都比较执行前后的 BRepCheck 错误并拒绝新增错误。
        bool validateAfterEachStage = true;

        // 为 true 时，所有阶段结束后再与整条流水线的原始输入比较一次。
        bool validateFinalResult = true;
    };

    // 单个模块执行器返回的统一结果。
    struct StageExecutionResult
    {
        TopoDS_Shape resultShape;
        bool modified = false;
        int repairCount = 0;
        std::vector<std::string> messages;
    };

    using StageExecutor = std::function<StageExecutionResult(const TopoDS_Shape&)>;

    // 一个可加入流水线的修复阶段；添加顺序就是执行顺序。
    struct Stage
    {
        std::string name;
        bool enabled = true;
        StageExecutor execute;
    };

    // 单个阶段的执行和提交情况，供 ApplicationInit 输出统一日志。
    struct StageReport
    {
        std::string name;
        bool enabled = true;
        bool executed = false;
        bool accepted = false;
        bool modified = false;
        int repairCount = 0;
        TopAbs_ShapeEnum inputShapeType = TopAbs_SHAPE;
        TopAbs_ShapeEnum outputShapeType = TopAbs_SHAPE;
        std::string rejectionReason;
        std::vector<std::string> messages;
    };

    // 一次完整流水线运行的结果。
    struct Result
    {
        TopoDS_Shape resultShape;
        bool modified = false;
        bool rolledBack = false;
        int totalRepairCount = 0;
        std::vector<StageReport> stageReports;
        std::vector<std::string> messages;
    };

    explicit ShapeRepairPipeline(const Params& params = Params());

    // 添加自定义修复阶段；执行器不得在失败或 modified=false 时原地破坏输入 shape。
    // 后续模块只需提供符合 StageExecutor 的适配函数，无需修改流水线核心。
    void addStage(const Stage& stage);
    void addStage(Stage&& stage);

    void clearStages();
    bool empty() const { return stages_.empty(); }

    const std::vector<Stage>& getStages() const { return stages_; }

    // 对一个独立拓扑容器运行所有已启用阶段。
    Result perform(const TopoDS_Shape& inputShape) const;

private:
    Params params_;
    std::vector<Stage> stages_;
};
