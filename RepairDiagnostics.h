#pragma once

#include <string>
#include <vector>

#include <BRepCheck_Status.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>

class BRepCheck_Analyzer;

// 统一封装修复算法使用的 BRepCheck 差异分析和结构化日志格式。
class RepairDiagnostics
{
public:
    struct BRepCheckFailureInfo
    {
        // 发生错误的拓扑子形状，用于比较修复前后是否为同一对象。
        TopoDS_Shape shape;

        // 子形状类型，用于避免不同拓扑层级之间误匹配。
        TopAbs_ShapeEnum type = TopAbs_SHAPE;

        // 同一子形状上所有非 NoError 的检查状态。
        std::vector<BRepCheck_Status> statuses;
    };

    using FailureList = std::vector<BRepCheckFailureInfo>;

    // 将 BRepCheck 状态转换为稳定、可读的英文名称。
    static const char* BRepCheckStatusName(BRepCheck_Status status);

    // 收集 shape 中所有 BRepCheck 错误；异常信息通过可选参数返回。
    static FailureList CollectBRepCheckFailures(
        const TopoDS_Shape& shape,
        std::string* errorMessage = nullptr);

    // 判断一个检查错误是否已经存在于修复前的错误列表中。
    static bool IsPreExistingBRepCheckFailure(
        const BRepCheckFailureInfo& failure,
        const FailureList& existingFailures);

    // 返回修复后新出现、修复前不存在的 BRepCheck 错误。
    static FailureList FindNewBRepCheckFailures(
        const FailureList& existingFailures,
        const FailureList& resultFailures);

    // 按 shape 类型和状态的出现次数比较错误，适用于修复后 TShape 身份发生变化的场景。
    static FailureList FindBRepCheckFailureCountIncrease(
        const FailureList& existingFailures,
        const FailureList& resultFailures);

    // 输出已收集错误的简洁说明。
    static std::string DescribeBRepCheckFailures(
        const FailureList& failures,
        int maxReported = 30);

    // 输出指定 analyzer 中带临时子形状序号的错误说明。
    static std::string DescribeBRepCheckFailures(
        BRepCheck_Analyzer& analyzer,
        const TopoDS_Shape& shape,
        int maxReported = 30);

    // 将连续的失败原因拆分成更适合日志阅读的多行文本。
    static std::string MakeFailureReasonReadable(const std::string& reason);

    // 为多行文本的每一行添加统一缩进。
    static std::string IndentMultiline(
        const std::string& text,
        const std::string& indent);

    // 生成统一的策略失败日志块。
    static std::string StrategyFailureBlock(
        const std::string& strategyName,
        const std::string& reason);

    // 添加带可选空行的迭代标题。
    static void AppendIterationHeader(
        std::vector<std::string>& messages,
        int iteration,
        bool addLeadingBlankLine);

    // 添加输入或结果 shape 的整体 BRepCheck 摘要。
    static void AppendShapeCheckSummary(
        std::vector<std::string>& messages,
        const std::string& shapeLabel,
        const TopoDS_Shape& shape);
};
