#include "Logger.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <initializer_list>
#include <sstream>
#include <system_error>
#include <utility>

namespace
{
    bool startsWith(const std::string& value, const std::string& prefix)
    {
        return value.size() >= prefix.size()
            && value.compare(0, prefix.size(), prefix) == 0;
    }

    bool containsAny(
        const std::string& value,
        std::initializer_list<const char*> needles)
    {
        for (const char* needle : needles) {
            if (value.find(needle) != std::string::npos) {
                return true;
            }
        }
        return false;
    }
}

Logger::Logger(
    const std::filesystem::path& filePath,
    std::string module,
    Level minimumLevel)
    : filePath_(filePath),
      module_(std::move(module)),
      minimumLevel_(minimumLevel)
{
    // 日志只是诊断辅助功能。即使目录创建失败，也不能抛出异常或中断几何修复；
    // 调用方可以通过 isOpen() 检查结果，并把日志故障输出到标准错误流。
    std::error_code directoryError;
    const std::filesystem::path parent = filePath_.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent, directoryError);
    }

    if (!directoryError) {
        stream_.open(filePath_, std::ios::out | std::ios::trunc);
    }
}

bool Logger::isOpen() const
{
    return stream_.is_open();
}

bool Logger::shouldLog(Level level) const
{
    // 等级值按 Debug < Info < Warn < Error < Off 排列，因此可以直接比较阈值。
    return minimumLevel_ != Level::Off
        && level >= minimumLevel_
        && level < Level::Off;
}

void Logger::run(
    Level level,
    const std::string& phase,
    const std::string& fields)
{
    write(level, "run", phase, fields);
}

void Logger::container(
    Level level,
    std::size_t containerIndex,
    const std::string& phase,
    const std::string& fields)
{
    write(
        level,
        "container",
        phase,
        appendField(integerField("container", static_cast<long long>(containerIndex)), fields));
}

void Logger::containerText(
    Level level,
    std::size_t containerIndex,
    const std::string& phase,
    const std::string& fields,
    const std::string& textName,
    const std::string& text)
{
    const std::string containerFields = appendField(
        integerField("container", static_cast<long long>(containerIndex)),
        fields);

    if (containsLineBreak(text)) {
        writeText(level, "container", phase, containerFields, textName, text);
        return;
    }

    write(level, "container", phase, appendField(containerFields, textField(textName, text)));
}

void Logger::containerMessages(
    std::size_t containerIndex,
    const std::string& algorithm,
    const std::vector<std::string>& messages)
{
    if (messages.empty()) {
        // “算法没有返回诊断内容”只对排查问题有帮助，默认 INFO 模式不记录，
        // 避免为正常 Container 增加无意义的日志行。
        container(
            Level::Debug,
            containerIndex,
            "process",
            appendField(
                appendField(textField("stage", "diagnostic"), textField("action", "no_detail")),
                textField("algorithm", algorithm)));
        return;
    }

    // 现有算法返回的是非结构化文字。统一在这一处完成等级和阶段归类，
    // 可以继续兼容旧消息接口，同时避免每个算法各自增加文件写入器和模式开关。
    for (const std::string& message : messages) {
        if (message.empty()) {
            continue;
        }

        std::string fields = textField("stage", classifyMessageStage(message));
        fields = appendField(fields, textField("action", classifyMessageAction(message)));
        fields = appendField(fields, textField("algorithm", algorithm));
        containerText(
            classifyMessageLevel(message),
            containerIndex,
            "process",
            fields,
            "message",
            message);
    }
}

Logger::Level Logger::levelFromEnvironment(Level fallback)
{
    // 构造 Logger 时只调用一次。运行中修改环境变量只会影响下一个 Logger，
    // 不会改变当前正在写入的任务。
    std::string configuredLevel;
#ifdef _WIN32
    char* rawValue = nullptr;
    std::size_t valueLength = 0;
    if (_dupenv_s(&rawValue, &valueLength, "OCC_HEAL_LOG_LEVEL") == 0
        && rawValue != nullptr) {
        configuredLevel.assign(rawValue);
    }
    std::free(rawValue);
#else
    const char* rawValue = std::getenv("OCC_HEAL_LOG_LEVEL");
    if (rawValue != nullptr) {
        configuredLevel.assign(rawValue);
    }
#endif
    return configuredLevel.empty() ? fallback : parseLevel(configuredLevel, fallback);
}

Logger::Level Logger::parseLevel(const std::string& value, Level fallback)
{
    const std::string normalized = lower(value);
    if (normalized == "debug") {
        return Level::Debug;
    }
    if (normalized == "info") {
        return Level::Info;
    }
    if (normalized == "warn" || normalized == "warning") {
        return Level::Warn;
    }
    if (normalized == "error") {
        return Level::Error;
    }
    if (normalized == "off" || normalized == "none") {
        return Level::Off;
    }
    return fallback;
}

const char* Logger::levelName(Level level)
{
    switch (level) {
    case Level::Debug:
        return "DEBUG";
    case Level::Info:
        return "INFO";
    case Level::Warn:
        return "WARN";
    case Level::Error:
        return "ERROR";
    case Level::Off:
        return "OFF";
    }
    return "INFO";
}

std::string Logger::textField(const std::string& name, const std::string& value)
{
    return name + "=" + quote(value);
}

std::string Logger::pathField(
    const std::string& name,
    const std::filesystem::path& value)
{
    // generic_string() 固定使用“/”；lexically_normal() 同时清理多余的
    // “.”、“..”和重复分隔符，但不会访问文件系统或要求目标路径已经存在。
    return textField(name, value.lexically_normal().generic_string());
}

std::string Logger::boolField(const std::string& name, bool value)
{
    return name + "=" + (value ? "true" : "false");
}

std::string Logger::integerField(const std::string& name, long long value)
{
    return name + "=" + std::to_string(value);
}

std::string Logger::numberField(const std::string& name, double value)
{
    std::ostringstream stream;
    stream << name << '=' << std::setprecision(12) << value;
    return stream.str();
}

std::string Logger::appendField(const std::string& fields, const std::string& field)
{
    if (fields.empty()) {
        return field;
    }
    if (field.empty()) {
        return fields;
    }
    return fields + " " + field;
}

void Logger::write(
    Level level,
    const char* eventName,
    const std::string& phase,
    const std::string& fields)
{
    if (!shouldLog(level) || !stream_.is_open()) {
        return;
    }

    // 在同一个锁内写完一整条记录，防止并发调用把不同日志的前缀和字段交叉写入。
    std::lock_guard<std::mutex> lock(mutex_);
    writePrefixUnlocked(level, eventName, phase, fields);
    stream_ << '\n';
    // 第三方几何计算可能在执行过程中异常退出。每写完一条有效日志就立即刷新，
    // 可以尽量保留故障发生前最后一段诊断信息。
    stream_.flush();
}

void Logger::writeText(
    Level level,
    const char* eventName,
    const std::string& phase,
    const std::string& fields,
    const std::string& textName,
    const std::string& text)
{
    if (!shouldLog(level) || !stream_.is_open()) {
        return;
    }

    // 首行和所有正文行必须在同一临界区内写完，否则不同线程的日志块可能交叉。
    std::lock_guard<std::mutex> lock(mutex_);
    writePrefixUnlocked(level, eventName, phase, fields);
    stream_ << ' ' << textName << ":\n";
    writeIndentedTextUnlocked(text);
    stream_.flush();
}

void Logger::writePrefixUnlocked(
    Level level,
    const char* eventName,
    const std::string& phase,
    const std::string& fields)
{
    stream_ << '[' << levelName(level) << "] [" << timestampNow() << "] ["
            << module_ << "] event=" << eventName
            << " phase=" << phase;
    if (!fields.empty()) {
        stream_ << ' ' << fields;
    }
}

void Logger::writeIndentedTextUnlocked(const std::string& text)
{
    // 每行统一增加四个空格，同时兼容 Windows(CRLF)、Linux(LF) 和旧式 CR 换行。
    // 保留算法原本的行内缩进，使嵌套失败原因在普通文本编辑器中仍有层次。
    std::size_t lineStart = 0;
    while (lineStart < text.size()) {
        std::size_t lineEnd = lineStart;
        while (lineEnd < text.size() && text[lineEnd] != '\r' && text[lineEnd] != '\n') {
            ++lineEnd;
        }

        stream_ << "    " << text.substr(lineStart, lineEnd - lineStart) << '\n';
        if (lineEnd == text.size()) {
            break;
        }

        if (text[lineEnd] == '\r'
            && lineEnd + 1 < text.size()
            && text[lineEnd + 1] == '\n') {
            lineStart = lineEnd + 2;
        }
        else {
            lineStart = lineEnd + 1;
        }
    }
}

Logger::Level Logger::classifyMessageLevel(const std::string& message)
{
    const std::string normalized = lower(message);

    // 判断顺序不能随意调整：最终失败和可恢复异常必须优先于基于关键字的
    // DEBUG 归类，否则一条包含内部计算词汇的错误消息可能被错误降级。
    if (containsAny(
            normalized,
            {"input shape is null", "result shape is null", "all repair strategies failed"})) {
        return Level::Error;
    }

    if (startsWith(normalized, "warning:")
        || containsAny(normalized, {"rejected", "could not", "cannot "})) {
        return Level::Warn;
    }

    if (containsAny(
            normalized,
            {
                "parameters:",
                "candidate skipped:",
                "stage disabled by pipeline configuration",
                "correspondence node ",
                "projection correspondence midline topology",
                "chaina edge ",
                "chainb edge ",
                "wire traversal check",
                "adjacency after",
                "adjacency before",
                "sample count",
                "sample_count",
                "mapping",
                "directedgeuses="
            })) {
        return Level::Debug;
    }

    if (containsAny(normalized, {"has existing brepcheck failures", "has brepcheck failures"})) {
        return Level::Warn;
    }

    return Level::Info;
}

const char* Logger::classifyMessageStage(const std::string& message)
{
    const std::string normalized = lower(message);
    if (containsAny(normalized, {"parameters:", "iteration limit"})) {
        return "config";
    }
    if (containsAny(normalized, {"brepcheck", "shape has", "shape check", "adjacency"})) {
        return "validate";
    }
    if (containsAny(normalized, {"iteration ", "candidate count", "candidate skipped"})) {
        return "detect";
    }
    return "repair";
}

const char* Logger::classifyMessageAction(const std::string& message)
{
    const std::string normalized = lower(message);
    if (containsAny(normalized, {"parameters:", "iteration limit"})) {
        return "configuration";
    }
    if (containsAny(normalized, {"brepcheck", "shape has", "shape check"})) {
        return "shape_check";
    }
    if (normalized.find("candidate count") != std::string::npos) {
        return "candidate_detection";
    }
    if (normalized.find("iteration ") != std::string::npos) {
        return "iteration";
    }
    if (containsAny(normalized, {"repaired candidate", "repair accepted", "repair succeeded"})) {
        return "succeeded";
    }
    if (normalized.find("repaired by") != std::string::npos
        && normalized.find("earlier strateg") != std::string::npos) {
        return "fallback_succeeded";
    }
    if (containsAny(normalized, {"failed", "rejected"})) {
        return "failed";
    }
    return "diagnostic";
}

std::string Logger::timestampNow()
{
    const auto now = std::chrono::system_clock::now();
    const auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    const std::time_t currentTime = std::chrono::system_clock::to_time_t(now);

    std::tm localTime{};
#ifdef _WIN32
    localtime_s(&localTime, &currentTime);
#else
    localtime_r(&currentTime, &localTime);
#endif

    std::ostringstream stream;
    stream << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S")
           << '.' << std::setw(3) << std::setfill('0') << milliseconds.count();
    return stream.str();
}

bool Logger::containsLineBreak(const std::string& value)
{
    return value.find_first_of("\r\n") != std::string::npos;
}

std::string Logger::quote(const std::string& value)
{
    // 普通结构化字段仍要转义引号、反斜杠和控制字符。需要展示的多行正文应通过
    // containerText() 写入；这里保留换行转义作为误用普通字段时的安全兜底。
    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');
    for (char ch : value) {
        switch (ch) {
        case '\\':
            escaped += "\\\\";
            break;
        case '"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped.push_back(ch);
            break;
        }
    }
    escaped.push_back('"');
    return escaped;
}

std::string Logger::lower(const std::string& value)
{
    std::string normalized = value;
    std::transform(
        normalized.begin(),
        normalized.end(),
        normalized.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return normalized;
}
