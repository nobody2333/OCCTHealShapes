#pragma once

#include <string>
#include <vector>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>

// SewShells 用于把同一拓扑父容器 compound 下的多个 shell 做 sewing 合并。
// 该类只负责 OCCT 拓扑形状层面的 sewing，不处理 XCAF label、颜色、名称等文档属性。
class SewShells {
public:
    // sewing 参数。默认策略偏保守：只接受没有新增明显拓扑问题的缝合结果。
    struct Params {
        // 流水线保守预处理：按输入 shell 归属局部缝合，不丢面、不处理闭壳。
        bool conservative = false;
        // BRepBuilderAPI_Sewing 使用的连接容差。
        double sewingTolerance = 1.0e-4;
        // sewing 允许使用的最小局部容差。
        double minTolerance = 1.0e-7;
        // sewing 允许使用的最大局部容差，避免容差被异常放大。
        double maxTolerance = 1.0e-2;
        // 同一父 compound 下至少有多少个直接 shell 子节点时才执行 sewing。
        int minShellCountPerGroup = 2;
        // 是否让 sewing 更新 same-parameter 状态。
        bool sameParameterMode = true;
        // 是否启用 face sewing。
        bool faceMode = true;
        // 是否处理 floating edge (不属于面的游离边）
        bool floatingEdgesMode = false;
        // 是否考虑局部 tolerance。
        bool localTolerancesMode = true;
        // 是否允许 non-manifold sewing。
        bool nonManifoldMode = false;
        // sewing 后是否对生成的 shell 做 face orientation 修复。
        bool fixShellOrientation = true;
        // 为 true 时要求一组 shell 必须缝合为单个 shell；为 false 时 shell 数量减少即可接受。
        bool requireSingleShell = false;
        // 为 true 时，在 sewing 后使用 ShapeFix_Shape 修复 sewing 引入的整体拓扑问题。
        bool runShapeFixAfterSewing = false;
        // 为 true 时允许丢弃 sewing 结果顶层裸露的 FACE/WIRE/EDGE/VERTEX；默认拒绝这类结果，避免静默丢面。
        bool discardLooseTopLevelSubshapesAfterSewing = false;
        // 如果 sewing 后 free edge 数量增加，是否拒绝本次结果。
        bool rejectFreeEdgeIncrease = true;
        // 如果 sewing 后 BRepCheck 不通过，是否拒绝本次结果。
        bool rejectInvalidResult = true;
        // 调试开关；当前关键过程始终写入 messages_。
        bool verbose = false;
    };

    explicit SewShells(const Params& params = Params());

    // tolerance 使用模型单位，独立于细缝宽度参数。
    static Params ConservativeParams(double tolerance = 1.0e-4);

    // 设置输入 shape，并清空上一轮执行状态。
    void setInputShape(const TopoDS_Shape& shape);

    // 获取最近一次执行后的结果；没有成功 sewing 时等于输入 shape。
    TopoDS_Shape getResultShape() const { return resultShape_; }

    // 获取最近一次设置的输入 shape。
    TopoDS_Shape getInputShape() const { return inputShape_; }

    // 最近一次执行是否发生了实际 sewing。
    bool isModified() const { return modified_; }

    // 最近一次执行成功缝合的 shell group 数量。
    int getSewedGroupCount() const { return sewedGroupCount_; }

    // 最近一次执行过程日志。
    const std::vector<std::string>& getMessages() const { return messages_; }

    // 使用当前 inputShape_ 执行 sewing，返回成功缝合的 shell group 数量。
    int perform();

    // 便捷入口：先 setInputShape(shape)，再 perform()。
    int perform(const TopoDS_Shape& shape);

private:
    // 对 shape 递归处理；只在 compound 的直接 shell 子节点之间执行 sewing。
    bool sewShapeRecursively(
        const TopoDS_Shape& sourceShape,
        TopoDS_Shape& resultShape,
        int& sewedGroupCount);

    // 对同一父 compound 下收集到的一组 shell 执行 sewing。
    bool sewShellGroup(
        const std::vector<TopoDS_Shell>& shellGroup,
        TopoDS_Shape& sewedShape,
        std::string& errorMessage);

    // 保守模式：固定输入 shell 归属，仅试缝跨壳自由边所在的两个面。
    bool sewCrossShellBoundaries(
        const std::vector<TopoDS_Shell>& shellGroup,
        TopoDS_Shape& sewedShape,
        std::string& errorMessage);

    // 构造只包含当前 shell group 的临时 compound，用于统计和校验。
    TopoDS_Shape makeShellGroupCompound(const std::vector<TopoDS_Shell>& shellGroup) const;

    // 对 sewing 结果中的 shell 做局部方向修复；compound 会递归处理其直接子 shape。
    TopoDS_Shape fixShellOrientationsInShape(
        const TopoDS_Shape& shape,
        std::string& warningMessage) const;

    // 根据 requireSingleShell 和 shell 数量变化决定是否接受 sewing 结果。
    bool acceptSewedResult(
        const TopoDS_Shape& rawSewedShape,
        int inputShellCount,
        TopoDS_Shape& acceptedShape,
        int& outputShellCount,
        std::string& errorMessage) const;

    // 对 sewing 后的 shape 应用 ShapeFix_Shape；用于修复 UnorientableShape 等整体拓扑问题。
    bool applyShapeFixAfterSewing(
        TopoDS_Shape& shape,
        const std::string& context,
        std::string& errorMessage);

    // 检查并按参数处理 sewing 结果顶层裸露的 FACE/WIRE/EDGE/VERTEX。
    bool handleLooseTopLevelSubshapes(
        TopoDS_Shape& shape,
        const std::string& context,
        std::string& errorMessage);

    // 判断 shape 是否是不能作为 shell sewing 结果顶层直接保留的裸拓扑。
    bool isLooseTopLevelSubshape(const TopoDS_Shape& shape) const;

    // 判断 shape 顶层是否包含裸拓扑；只检查直接子节点，不递归进入 compound。
    bool hasLooseTopLevelSubshape(const TopoDS_Shape& shape) const;

    // 从 sewing 结果中提取唯一 shell；严格模式下用于确认结果已经缝合成单个 shell。
    bool extractSingleShell(
        const TopoDS_Shape& shape,
        TopoDS_Shell& shell,
        std::string& errorMessage) const;

    // 校验 sewing 结果没有明显变差。
    bool validateSewingResult(
        const TopoDS_Shape& beforeShape,
        const TopoDS_Shape& afterShape,
        std::string& errorMessage,
        std::string& warningMessage) const;

    // 统计 shape 内 shell 数量。
    int countShells(const TopoDS_Shape& shape) const;

    // 统计 shape 的直接子拓扑数量；非 compound 按 1 个直接组件处理。
    int countDirectChildren(const TopoDS_Shape& shape) const;

    // 输出 shape 直接子拓扑类型分布，用于判断 sewing 是否只是产生了包装层。
    std::string describeDirectChildren(const TopoDS_Shape& shape) const;

    // 追加过程日志。
    void appendMessage(const std::string& message);

private:
    // 当前 sewing 参数。
    Params params_;

    // 当前输入 shape。
    TopoDS_Shape inputShape_;

    // 最近一次执行后的结果 shape。
    TopoDS_Shape resultShape_;

    // 最近一次执行是否修改了 shape。
    bool modified_ = false;

    // 最近一次成功 sewing 的 shell group 数量。
    int sewedGroupCount_ = 0;

    // 最近一次执行过程日志。
    std::vector<std::string> messages_;
};
