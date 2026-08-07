#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <vector>

// 所有几何修复入口共用的日志写入器。
//
// 一级事件只保留两种，避免事件名称随着算法增多而不断膨胀：
//   event=run       - 一次完整的修复和导出任务
//   event=container - 该任务处理的一个可编辑 solid/shell 容器
//
// 更具体的含义由 phase/stage/action/status 字段表达。除非调用方显式指定，
// 每个 Logger 创建时只读取一次 OCC_HEAL_LOG_LEVEL，作为本次任务的最低输出等级。
class Logger
{
public:
    // 枚举值必须按照“信息最详细”到“问题最严重”的顺序排列，
    // shouldLog() 会直接比较枚举值来判断一条日志是否达到输出阈值。
    enum class Level
    {
        Debug = 0,
        Info = 1,
        Warn = 2,
        Error = 3,
        Off = 4
    };

    static Level levelFromEnvironment(Level fallback = Level::Info);
    static Level parseLevel(const std::string& value, Level fallback = Level::Info);
    static const char* levelName(Level level);

    Logger(
        const std::filesystem::path& filePath,
        std::string module,
        Level minimumLevel = levelFromEnvironment());

    bool isOpen() const;
    const std::filesystem::path& filePath() const { return filePath_; }
    Level minimumLevel() const { return minimumLevel_; }
    // 应在任务开始前设置最低输出等级。一个 Logger 通常只服务于一个修复入口，
    // 不建议在任务执行过程中动态修改等级，以免同一份日志的详细程度前后不一致。
    void setMinimumLevel(Level level) { minimumLevel_ = level; }
    bool shouldLog(Level level) const;

    void run(Level level, const std::string& phase, const std::string& fields = {});
    void container(
        Level level,
        std::size_t containerIndex,
        const std::string& phase,
        const std::string& fields = {});

    // 写入带文本正文的 Container 日志。短文本仍作为普通字段写在首行；
    // 多行文本固定使用“字段名:”加缩进正文的易读格式，不提供格式模式开关。
    void containerText(
        Level level,
        std::size_t containerIndex,
        const std::string& phase,
        const std::string& fields,
        const std::string& textName,
        const std::string& text);

    // 兼容当前算法仍然通过 vector<string> 返回诊断内容的接口。
    // perform() 结束后，本方法才会逐条判断等级并写入时间戳，因此这里的时间表示
    // “消息写入日志的时间”，不一定是算法内部计算实际发生的精确时刻。
    void containerMessages(
        std::size_t containerIndex,
        const std::string& algorithm,
        const std::vector<std::string>& messages);

    // 字段辅助函数负责为单行值统一添加引号和转义字符。确定是单行的路径、状态、
    // 错误原因等应通过这些函数生成；可能包含换行的正文应使用 containerText()。
    static std::string textField(const std::string& name, const std::string& value);
    // 路径字段统一使用“/”作为分隔符，避免 Windows 路径在日志中混用“/”和“\”。
    static std::string pathField(
        const std::string& name,
        const std::filesystem::path& value);
    static std::string boolField(const std::string& name, bool value);
    static std::string integerField(const std::string& name, long long value);
    static std::string numberField(const std::string& name, double value);
    static std::string appendField(const std::string& fields, const std::string& field);

private:
    void write(
        Level level,
        const char* eventName,
        const std::string& phase,
        const std::string& fields);
    void writeText(
        Level level,
        const char* eventName,
        const std::string& phase,
        const std::string& fields,
        const std::string& textName,
        const std::string& text);
    void writePrefixUnlocked(
        Level level,
        const char* eventName,
        const std::string& phase,
        const std::string& fields);
    void writeIndentedTextUnlocked(const std::string& text);

    static Level classifyMessageLevel(const std::string& message);
    static const char* classifyMessageStage(const std::string& message);
    static const char* classifyMessageAction(const std::string& message);
    static std::string timestampNow();
    static bool containsLineBreak(const std::string& value);
    static std::string quote(const std::string& value);
    static std::string lower(const std::string& value);

    std::filesystem::path filePath_;
    std::string module_;
    Level minimumLevel_ = Level::Info;
    mutable std::mutex mutex_;
    std::ofstream stream_;
};
