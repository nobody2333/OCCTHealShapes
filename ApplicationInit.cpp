#include "ApplicationInit.h"

#include <cctype>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <IGESCAFControl_Reader.hxx>
#include <IGESData_BasicEditor.hxx>
#include <IGESData_GlobalSection.hxx>
#include <IGESData_IGESModel.hxx>
#include <Message.hxx>
#include <Message_Printer.hxx>
#include <Message_SequenceOfPrinters.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDF_LabelSequence.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TColStd_SequenceOfAsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

#include "Logger.h"
#include "TopologyUtils.h"

namespace {
    using RepairClock = std::chrono::steady_clock;

    // Logger 接收一段扁平的字段字符串。统一通过 Logger 的字段辅助函数构造，
    // 确保所有文本都使用相同的引号和转义规则。
    std::string logFields(std::initializer_list<std::string> fields)
    {
        std::string result;
        for (const std::string& field : fields) {
            result = Logger::appendField(result, field);
        }
        return result;
    }

    long long elapsedMilliseconds(const RepairClock::time_point& start)
    {
        // steady_clock 不受系统时间校准影响，适合计算耗时较长的几何任务时长。
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            RepairClock::now() - start).count();
    }

    void reportLogOpenFailure(const Logger& logger)
    {
        // 日志打开失败不能影响修复任务本身；几何修复继续执行，同时在控制台
        // 明确报告日志文件不可用。
        if (!logger.isOpen()) {
            std::cerr << "Warning: Cannot open repair log file: "
                      << logger.filePath().string() << std::endl;
        }
    }

    std::string shapeTypeField(const TopoDS_Shape& shape)
    {
        return Logger::textField("shape_type", TopologyUtils::ShapeTypeName(shape));
    }

    void logRunFailure(
        Logger& logger,
        const RepairClock::time_point& start,
        const std::string& stage,
        const std::string& reason)
    {
        // 所有提前返回路径都通过同一函数写 Run 失败汇总，避免日志只有
        // event=run phase=start，却没有与之对应的 phase=end。
        logger.run(
            Logger::Level::Error,
            "end",
            logFields({
                Logger::textField("status", "failed"),
                Logger::textField("stage", stage),
                Logger::textField("reason", reason),
                Logger::integerField("duration_ms", elapsedMilliseconds(start))
            }));
    }

    std::string labelEntryString(const TDF_Label& label)
    {
        if (label.IsNull()) {
            return "<null>";
        }

        std::ostringstream stream;
        label.EntryDump(stream);
        return stream.str();
    }

    struct ShapeLabelDebugInfo {
        TDF_Label label;
        TDF_Label parentLabel;
        TopAbs_ShapeEnum shapeType;
    };

    void collectShellCompoundLabels(const TDF_Label& label,
                                    std::vector<ShapeLabelDebugInfo>& shellLabels,
                                    std::vector<ShapeLabelDebugInfo>& compoundLabels)
    {
        if (XCAFDoc_ShapeTool::IsShape(label)) {
            const TopoDS_Shape shape = XCAFDoc_ShapeTool::GetShape(label);
            if (!shape.IsNull()) {
                const TopAbs_ShapeEnum shapeType = shape.ShapeType();
                if (shapeType == TopAbs_SHELL) {
                    shellLabels.push_back({ label, label.Father(), shapeType });
                }
                else if (shapeType == TopAbs_COMPOUND) {
                    compoundLabels.push_back({ label, label.Father(), shapeType });
                }
            }
        }

        for (TDF_ChildIterator childIterator(label); childIterator.More(); childIterator.Next()) {
            collectShellCompoundLabels(childIterator.Value(), shellLabels, compoundLabels);
        }
    }

    void appendUniqueParentLabel(const TDF_Label& parentLabel, std::vector<TDF_Label>& parentLabels)
    {
        for (const TDF_Label& existingParentLabel : parentLabels) {
            if (existingParentLabel == parentLabel) {
                return;
            }
        }

        parentLabels.push_back(parentLabel);
    }

    void dumpShellCompoundParentDebug(const TDF_Label& rootLabel)
    {
        std::vector<ShapeLabelDebugInfo> shellLabels;
        std::vector<ShapeLabelDebugInfo> compoundLabels;
        collectShellCompoundLabels(rootLabel, shellLabels, compoundLabels);

        std::vector<TDF_Label> parentLabels;
        for (const ShapeLabelDebugInfo& shellInfo : shellLabels) {
            appendUniqueParentLabel(shellInfo.parentLabel, parentLabels);
        }
        for (const ShapeLabelDebugInfo& compoundInfo : compoundLabels) {
            appendUniqueParentLabel(compoundInfo.parentLabel, parentLabels);
        }

        std::cout << "[collectTopLevelTargets debug] shell count: " << shellLabels.size()
                  << ", compound count: " << compoundLabels.size()
                  << ", parent group count: " << parentLabels.size() << std::endl;

        for (const TDF_Label& parentLabel : parentLabels) {
            std::vector<std::string> shellEntries;
            std::vector<std::string> compoundEntries;

            for (const ShapeLabelDebugInfo& shellInfo : shellLabels) {
                if (shellInfo.parentLabel == parentLabel) {
                    shellEntries.push_back(labelEntryString(shellInfo.label));
                }
            }

            for (const ShapeLabelDebugInfo& compoundInfo : compoundLabels) {
                if (compoundInfo.parentLabel == parentLabel) {
                    compoundEntries.push_back(labelEntryString(compoundInfo.label));
                }
            }

            std::cout << "[collectTopLevelTargets debug] parent="
                      << labelEntryString(parentLabel)
                      << ", shells=" << shellEntries.size()
                      << ", compounds=" << compoundEntries.size() << std::endl;

            for (const std::string& shellEntry : shellEntries) {
                std::cout << "  shell label: " << shellEntry << std::endl;
            }

            for (const std::string& compoundEntry : compoundEntries) {
                std::cout << "  compound label: " << compoundEntry << std::endl;
            }
        }
    }

    std::string normalizeStepLengthUnit(const TCollection_AsciiString& unitName)
    {
        std::string unit = unitName.ToCString();
        for (char& ch : unit) {
            ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        }

        if (unit == "METRE" || unit == "METER" || unit == "M") {
            return "M";
        }
        if (unit == "MILLIMETRE" || unit == "MILLIMETER" || unit == "MM") {
            return "MM";
        }
        if (unit == "CENTIMETRE" || unit == "CENTIMETER" || unit == "CM") {
            return "CM";
        }
        if (unit == "INCH" || unit == "IN") {
            return "INCH";
        }
        if (unit == "FOOT" || unit == "FT") {
            return "FT";
        }
        if (unit == "MICROMETRE" || unit == "MICROMETER" || unit == "UM") {
            return "UM";
        }
        if (unit == "KILOMETRE" || unit == "KILOMETER" || unit == "KM") {
            return "KM";
        }

        return "MM";
    }

    void printStepShapeType(const TopoDS_Shape inputShape) {
        int index = 0;
        for (TopoDS_Iterator it(inputShape); it.More(); it.Next()) {
            const TopoDS_Shape child = it.Value();
            if (child.ShapeType() == TopAbs_COMPOUND) {
                std::cout << "child[" << index << "] type = "
                    << TopologyUtils::ShapeTypeName(child.ShapeType())
                    << std::endl;
				printStepShapeType(child);
                continue;
            }
            std::cout << "child[" << index << "] type = "
                << TopologyUtils::ShapeTypeName(child.ShapeType())
                << std::endl;

            ++index;
        }
    }

// 控制Transfer减少输出信息
class ScopedOccInfoSilencer {
public:
    ScopedOccInfoSilencer()
    {
        const Handle(Message_Messenger)& messenger = Message::DefaultMessenger();
        for (Message_SequenceOfPrinters::Iterator it(messenger->Printers()); it.More(); it.Next()) {
            const Handle(Message_Printer)& printer = it.Value();
            printers_.push_back(printer);
            traceLevels_.push_back(printer->GetTraceLevel());
            printer->SetTraceLevel(Message_Warning);
        }
    }

    ~ScopedOccInfoSilencer()
    {
        for (std::size_t i = 0; i < printers_.size(); ++i) {
            if (!printers_[i].IsNull()) {
                printers_[i]->SetTraceLevel(traceLevels_[i]);
            }
        }
    }

private:
    std::vector<Handle(Message_Printer)> printers_;
    std::vector<Message_Gravity> traceLevels_;
};
}


ApplicationInit::ApplicationInit()
{
    anApp_ = new TDocStd_Application();
    anApp_->NewDocument("BinXCAF", aDoc_);
}

ApplicationInit::~ApplicationInit()
{
}

bool ApplicationInit::importStep(const std::string& filePath)
{
    STEPCAFControl_Reader aReader;
    std::cout << "Reading STEP file: " << filePath << "......" << std::endl;

    const IFSelect_ReturnStatus status = aReader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to read STEP file." << std::endl;
        return false;
    }

    TColStd_SequenceOfAsciiString lengthUnits;
    TColStd_SequenceOfAsciiString angleUnits;
    TColStd_SequenceOfAsciiString solidAngleUnits;
    aReader.ChangeReader().FileUnits(lengthUnits, angleUnits, solidAngleUnits);
    if (!lengthUnits.IsEmpty()) {
        inputStepUnit_ = normalizeStepLengthUnit(lengthUnits.First());
    }
    else {
        inputStepUnit_ = "MM";
    }
    std::cout << "Input STEP length unit: " << inputStepUnit_ << std::endl;

    std::cout << "STEP file read successfully." << std::endl;
    return transferAndCacheTopLevelShapes(aReader.Transfer(aDoc_), "STEP");
}

bool ApplicationInit::importStepAsShape(const std::string& filePath)
{
    STEPControl_Reader reader;
    importedStepShape_.Nullify();

    std::cout << "Reading STEP shape with STEPControl_Reader: " << filePath << "......" << std::endl;

    const IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: STEPControl_Reader failed to read STEP file." << std::endl;
        return false;
    }

    TColStd_SequenceOfAsciiString lengthUnits;
    TColStd_SequenceOfAsciiString angleUnits;
    TColStd_SequenceOfAsciiString solidAngleUnits;
    reader.FileUnits(lengthUnits, angleUnits, solidAngleUnits);
    if (!lengthUnits.IsEmpty()) {
        inputStepUnit_ = normalizeStepLengthUnit(lengthUnits.First());
    }
    else {
        inputStepUnit_ = "MM";
    }
    std::cout << "Input STEP length unit: " << inputStepUnit_ << std::endl;

    if (reader.NbRootsForTransfer() <= 0) {
        std::cerr << "Error: STEPControl_Reader found no transferable roots." << std::endl;
        return false;
    }

    const Standard_Integer transferredRootCount = reader.TransferRoots();
    if (transferredRootCount <= 0) {
        std::cerr << "Error: STEPControl_Reader failed to transfer roots." << std::endl;
        return false;
    }

    importedStepShape_ = reader.OneShape();

    // DEBUG: 临时信息
    printStepShapeType(importedStepShape_);

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: STEPControl_Reader produced a null shape." << std::endl;
        return false;
    }

    std::cout << "STEPControl_Reader transferred roots: " << transferredRootCount << std::endl;
    std::cout << "STEP shape read successfully." << std::endl;
    return true;
}

bool ApplicationInit::importIges(const std::string& filePath)
{
    IGESCAFControl_Reader aReader;

    std::cout << "Reading IGES file: " << filePath << "......" << std::endl;

    const IFSelect_ReturnStatus status = aReader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to read IGES file." << std::endl;
        return false;
    }

    const Handle(IGESData_IGESModel) igesModel = aReader.IGESModel();
    if (!igesModel.IsNull()) {
        const IGESData_GlobalSection& globalSection = igesModel->GlobalSection();
        const Handle(TCollection_HAsciiString) unitName = globalSection.UnitName();
        if (!unitName.IsNull() && unitName->Length() > 0) {
            inputStepUnit_ = normalizeStepLengthUnit(unitName->String());
        }
        else {
            const Standard_CString unitFlagName =
                IGESData_BasicEditor::UnitFlagName(globalSection.UnitFlag());
            inputStepUnit_ = (unitFlagName != nullptr && unitFlagName[0] != '\0')
                ? normalizeStepLengthUnit(TCollection_AsciiString(unitFlagName))
                : "MM";
        }
    }
    else {
        inputStepUnit_ = "MM";
    }
    std::cout << "Input IGES length unit: " << inputStepUnit_ << std::endl;

    std::cout << "IGES file read successfully." << std::endl;
    return transferAndCacheTopLevelShapes(aReader.Transfer(aDoc_), "IGES");
}

bool ApplicationInit::exportStep(const TopoDS_Shape& shape, const std::string& filePath) const
{
    if (shape.IsNull()) {
        std::cerr << "Error: The shape to export is null." << std::endl;
        return false;
    }

    Handle(TDocStd_Document) exportDoc;
    Handle(TDocStd_Application) exportApp = new TDocStd_Application();
    exportApp->NewDocument("BinXCAF", exportDoc);

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(exportDoc->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to create XCAF shape tool for STEP export." << std::endl;
        return false;
    }

    shapeTool->AddShape(shape, Standard_False, Standard_True);
    shapeTool->UpdateAssemblies();

    return writeXcafDocumentStep(exportDoc, filePath, false);
}

bool ApplicationInit::exportStepShape(const TopoDS_Shape& shape, const std::string& filePath) const
{
    if (shape.IsNull()) {
        std::cerr << "Error: The shape to export with STEPControl_Writer is null." << std::endl;
        return false;
    }

    const std::filesystem::path outputPath(filePath);
    if (outputPath.has_parent_path()) {
        std::filesystem::create_directories(outputPath.parent_path());
    }

    std::cout << "Writing STEP shape with STEPControl_Writer to: " << filePath << "......" << std::endl;
    setOutputStepUnitFromInput();

    STEPControl_Writer writer;

	// 利用RAII机制在Transfer过程中暂时提升OCCT消息级别，减少输出信息干扰
    ScopedOccInfoSilencer occInfoSilencer;
    if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone) {
        std::cerr << "Error: Failed to transfer shape to STEPControl_Writer." << std::endl;
        return false;
    }

    if (writer.Write(filePath.c_str()) == IFSelect_RetDone) {
        std::cout << "STEP shape written successfully." << std::endl;
        return true;
    }

    std::cerr << "Error: Failed to write STEP shape." << std::endl;
    return false;
}

int ApplicationInit::mergeDuplicateVerticesOnStepShapeAndExport(
    const std::string& filePath,
    const MergeDuplicateVertex::Params& params)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexStepShape" / "MergeDuplicateVertexStepShapeLog.txt";
    Logger logger(mergeLogPath, "MergeDuplicateVertex");
    reportLogOpenFailure(logger);
    logger.run(
        Logger::Level::Info,
        "start",
        logFields({
            Logger::textField("input_scope", "step_shape"),
            Logger::pathField("output_step", filePath),
            Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
        }));

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: Cannot merge duplicate vertices because imported STEP shape is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "imported STEP shape is null");
        return 0;
    }

    const RepairClock::time_point containerStart = RepairClock::now();
    logger.container(
        Logger::Level::Info,
        0,
        "start",
        logFields({shapeTypeField(importedStepShape_), Logger::textField("scope", "step_shape")}));

    MergeDuplicateVertex duplicateVertexMerger(params);
    duplicateVertexMerger.setInputShape(importedStepShape_);
    const int mergedCount = duplicateVertexMerger.perform();
    logger.containerMessages(0, "MergeDuplicateVertex", duplicateVertexMerger.getMessages());

    TopoDS_Shape processedShape = duplicateVertexMerger.getResultShape();
    if (processedShape.IsNull()) {
        std::cerr << "Error: MergeDuplicateVertex produced a null STEP shape result." << std::endl;
        logger.container(
            Logger::Level::Error,
            0,
            "end",
            logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
        logRunFailure(logger, runStart, "result_validate", "result shape is null");
        return 0;
    }

    if (processedShape.ShapeType() != importedStepShape_.ShapeType()) {
        std::cerr << "Error: MergeDuplicateVertex changed STEP shape type from "
                  << static_cast<int>(importedStepShape_.ShapeType()) << " to "
                  << static_cast<int>(processedShape.ShapeType()) << "." << std::endl;
        logger.container(
            Logger::Level::Error,
            0,
            "end",
            logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
        logRunFailure(logger, runStart, "result_validate", "result shape type changed");
        return 0;
    }

    importedStepShape_ = processedShape;

    logger.container(
        Logger::Level::Info,
        0,
        "end",
        logFields({
            Logger::textField(
                "status",
                duplicateVertexMerger.isModified() ? "success" : "no_change"),
            Logger::boolField("modified", duplicateVertexMerger.isModified()),
            Logger::integerField("repair_count", mergedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));

    if (!exportStepShape(importedStepShape_, filePath)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP shape");
        return -1;
    }

    logger.run(
        Logger::Level::Info,
        "end",
        logFields({
            Logger::textField("status", "success"),
            Logger::integerField("container_count", 1),
            Logger::integerField("modified_containers", duplicateVertexMerger.isModified() ? 1 : 0),
            Logger::integerField("repair_count", mergedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
        }));

    return mergedCount;
}

int ApplicationInit::sewShellsOnStepShapeAndExport(
    const std::string& filePath,
    const SewShells::Params& params)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path sewLogPath =
        std::filesystem::path("assets") / "fixedSewShellsStepShape" / "SewShellsStepShapeLog.txt";
    Logger logger(sewLogPath, "SewShells");
    reportLogOpenFailure(logger);
    logger.run(
        Logger::Level::Info,
        "start",
        logFields({
            Logger::textField("input_scope", "step_shape"),
            Logger::pathField("output_step", filePath),
            Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
        }));

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: Cannot sew shells because imported STEP shape is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "imported STEP shape is null");
        return 0;
    }

    const RepairClock::time_point containerStart = RepairClock::now();
    logger.container(
        Logger::Level::Info,
        0,
        "start",
        logFields({shapeTypeField(importedStepShape_), Logger::textField("scope", "step_shape")}));

    SewShells shellSewer(params);
    shellSewer.setInputShape(importedStepShape_);
    const int sewedGroupCount = shellSewer.perform();
    logger.containerMessages(0, "SewShells", shellSewer.getMessages());

    TopoDS_Shape processedShape = shellSewer.getResultShape();
    if (processedShape.IsNull()) {
        std::cerr << "Error: SewShells produced a null STEP shape result." << std::endl;
        logger.container(
            Logger::Level::Error,
            0,
            "end",
            logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
        logRunFailure(logger, runStart, "result_validate", "result shape is null");
        return 0;
    }

    if (processedShape.ShapeType() != importedStepShape_.ShapeType()) {
        std::cerr << "Error: SewShells changed STEP shape type from "
                  << TopologyUtils::ShapeTypeName(importedStepShape_) << " to "
                  << TopologyUtils::ShapeTypeName(processedShape) << "." << std::endl;
        logger.container(
            Logger::Level::Error,
            0,
            "end",
            logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
        logRunFailure(logger, runStart, "result_validate", "result shape type changed");
        return 0;
    }

    importedStepShape_ = processedShape;

    logger.container(
        Logger::Level::Info,
        0,
        "end",
        logFields({
            Logger::textField("status", shellSewer.isModified() ? "success" : "no_change"),
            Logger::boolField("modified", shellSewer.isModified()),
            Logger::integerField("repair_count", sewedGroupCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));

    if (!exportStepShape(importedStepShape_, filePath)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP shape");
        return -1;
    }

    logger.run(
        Logger::Level::Info,
        "end",
        logFields({
            Logger::textField("status", "success"),
            Logger::integerField("container_count", 1),
            Logger::integerField("modified_containers", shellSewer.isModified() ? 1 : 0),
            Logger::integerField("repair_count", sewedGroupCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
        }));

    return sewedGroupCount;
}

void ApplicationInit::collectTopLevelTargets(std::vector<EditableContainerTarget>& targets) const {
    if (aDoc_.IsNull()) {
        return;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        return;
    }
#ifdef _DEBUG
    dumpShellCompoundParentDebug(aDoc_->Main());
#endif
    TDF_LabelSequence topShapes;
    shapeTool->GetShapes(topShapes);
    
    for (Standard_Integer i = 1; i <= topShapes.Length(); ++i) {
        TDF_Label label = topShapes.Value(i);
        if (!XCAFDoc_ShapeTool::IsSimpleShape(label)) {
            continue;
        }
        const TopoDS_Shape shape = shapeTool->GetShape(label);
        if (!shape.IsNull() && shapeContainsSolidOrShell(shape)) {
            targets.push_back({ label, shape });
        }
    }
}

int ApplicationInit::mergeShortEdgesOnLowestContainersAndExport(const std::string& filePath,
                                                                bool useManualLengthThreshold,
                                                                double manualLengthThreshold,
                                                                double angleThresholdDeg)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path mergeLogPath = std::filesystem::path("assets") / "fixedShortEdge" / "MergeShortEdgeLog.txt";
    Logger logger(mergeLogPath, "MergeShortEdge");
    reportLogOpenFailure(logger);
    logger.run(
        Logger::Level::Info,
        "start",
        logFields({
            Logger::pathField("output_step", filePath),
            Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
        }));
    logger.run(
        Logger::Level::Info,
        "process",
        logFields({
            Logger::textField("stage", "config"),
            Logger::boolField("manual_length_threshold", useManualLengthThreshold),
            Logger::numberField("length_threshold", manualLengthThreshold),
            Logger::numberField("angle_threshold_deg", angleThresholdDeg)
        }));

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot merge short edges because the XCAF document is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    //collectEditableContainerTargets(targets);
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int mergedCount = 0;
    int replacedCount = 0;
    int skippedCount = 0;
    int failedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            ++skippedCount;
            logger.container(
                Logger::Level::Warn,
                i,
                "end",
                logFields({
                    Logger::textField("status", "skipped"),
                    Logger::textField("reason", "owner label is null"),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            continue;
        }

        // 重新读取 label 当前 shape，避免前面容器写回后导致装配派生的
        // compound 已刷新，而这里仍使用旧缓存。
        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            ++skippedCount;
            logger.container(
                Logger::Level::Warn,
                i,
                "end",
                logFields({
                    Logger::textField("status", "skipped"),
                    Logger::textField("reason", "input shape is null"),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            continue;
        }

        logger.container(Logger::Level::Info, i, "start", shapeTypeField(originalShape));

        MergeShortEdge shortEdgeMerger;
        if (useManualLengthThreshold) {
            shortEdgeMerger.setLengthThreshold(manualLengthThreshold);
        }
        shortEdgeMerger.setAngleThresholdDeg(angleThresholdDeg);
        shortEdgeMerger.setInputShape(originalShape);
        const int containerMergedCount = shortEdgeMerger.performMerge();
        logger.containerMessages(i, "MergeShortEdge", shortEdgeMerger.getMessages());

        TopoDS_Shape processedShape = shortEdgeMerger.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            ++failedCount;
            logger.container(
                Logger::Level::Error,
                i,
                "end",
                logFields({
                    Logger::textField("status", "failed"),
                    Logger::boolField("modified", false),
                    Logger::textField("action", "keep_original"),
                    Logger::textField("reason", "result shape is null"),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            ++failedCount;
            logger.container(
                Logger::Level::Error,
                i,
                "end",
                logFields({
                    Logger::textField("status", "failed"),
                    Logger::boolField("modified", false),
                    Logger::textField("action", "keep_original"),
                    Logger::textField("reason", "result shape type changed"),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            continue;
        }

        // 整个叶子容器只写回一次。这样会牺牲该局部装配 label 内部的
        // 细粒度装配语义，但能保留上层 STEP/XCAF 父级关系，避免 Cadfix
        // 将 shell 下的 face 识别成 parentless face。
        shapeTool->SetShape(target.ownerLabel, processedShape);
        mergedCount += containerMergedCount;
        ++replacedCount;
        logger.container(
            Logger::Level::Info,
            i,
            "end",
            logFields({
                Logger::textField("status", containerMergedCount > 0 ? "success" : "no_change"),
                Logger::boolField("modified", containerMergedCount > 0),
                Logger::integerField("repair_count", containerMergedCount),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after short-edge merge: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(
        failedCount > 0 ? Logger::Level::Warn : Logger::Level::Info,
        "end",
        logFields({
            Logger::textField("status", failedCount > 0 ? "partial_success" : "success"),
            Logger::integerField("container_count", static_cast<long long>(targets.size())),
            Logger::integerField("modified_containers", replacedCount),
            Logger::integerField("skipped_containers", skippedCount),
            Logger::integerField("failed_containers", failedCount),
            Logger::integerField("repair_count", mergedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
        }));

    return mergedCount;
}

int ApplicationInit::mergeDuplicateVerticesOnLowestContainersAndExport(
    const std::string& filePath,
    const MergeDuplicateVertex::Params& params)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertex" / "MergeDuplicateVertexLog.txt";
    Logger logger(mergeLogPath, "MergeDuplicateVertex");
    reportLogOpenFailure(logger);
    logger.run(
        Logger::Level::Info,
        "start",
        logFields({
            Logger::pathField("output_step", filePath),
            Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
        }));

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot merge duplicate vertices because the XCAF document is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int mergedCount = 0;
    int replacedCount = 0;
    int skippedCount = 0;
    int failedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "owner label is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "input shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        logger.container(Logger::Level::Info, i, "start", shapeTypeField(originalShape));

        MergeDuplicateVertex duplicateVertexMerger(params);
        duplicateVertexMerger.setInputShape(originalShape);
        const int containerMergedCount = duplicateVertexMerger.perform();
        logger.containerMessages(i, "MergeDuplicateVertex", duplicateVertexMerger.getMessages());

        TopoDS_Shape processedShape = duplicateVertexMerger.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"),
                Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (!duplicateVertexMerger.isModified()) {
            logger.container(Logger::Level::Info, i, "end", logFields({
                Logger::textField("status", "no_change"),
                Logger::boolField("modified", false),
                Logger::integerField("repair_count", 0),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        mergedCount += containerMergedCount;
        ++replacedCount;
        logger.container(Logger::Level::Info, i, "end", logFields({
            Logger::textField("status", "success"),
            Logger::boolField("modified", true),
            Logger::integerField("repair_count", containerMergedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after duplicate-vertex merge: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(
        failedCount > 0 ? Logger::Level::Warn : Logger::Level::Info,
        "end",
        logFields({
            Logger::textField("status", failedCount > 0 ? "partial_success" : "success"),
            Logger::integerField("container_count", static_cast<long long>(targets.size())),
            Logger::integerField("modified_containers", replacedCount),
            Logger::integerField("skipped_containers", skippedCount),
            Logger::integerField("failed_containers", failedCount),
            Logger::integerField("repair_count", mergedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
        }));

    return mergedCount;
}

int ApplicationInit::mergeDuplicateVerticesAndHealNarrowGapsOnLowestContainersAndExport(
    const std::string& filePath,
    const MergeDuplicateVertex::Params& duplicateVertexParams,
    const HealNarrowGap::Params& narrowGapParams)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path chainLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexNarrowGap" / "MergeDuplicateVertexHealNarrowGapLog.txt";
    Logger logger(chainLogPath, "MergeDuplicateVertex+HealNarrowGap");
    reportLogOpenFailure(logger);
    logger.run(Logger::Level::Info, "start", logFields({
        Logger::pathField("output_step", filePath),
        Logger::textField("log_level", Logger::levelName(logger.minimumLevel())),
        Logger::integerField("configured_stage_count", 2)
    }));

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot run duplicate-vertex/narrow-gap repair because the XCAF document is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int duplicateMergedCount = 0;
    int narrowGapRepairedCount = 0;
    int replacedCount = 0;
    int skippedCount = 0;
    int warningCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"), Logger::textField("reason", "owner label is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"), Logger::textField("reason", "input shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        logger.container(Logger::Level::Info, i, "start", shapeTypeField(originalShape));

        TopoDS_Shape workingShape = originalShape;
        bool containerModified = false;
        bool containerWarning = false;

        // 每个通过校验的阶段都会推进 workingShape。如果某阶段声称已修改，
        // 但结果为空或 Shape 类型改变，则保留上一份有效 Shape，让后续阶段继续执行。
        MergeDuplicateVertex duplicateVertexMerger(duplicateVertexParams);
        duplicateVertexMerger.setInputShape(workingShape);
        const int containerDuplicateMergedCount = duplicateVertexMerger.perform();
        logger.containerMessages(i, "MergeDuplicateVertex", duplicateVertexMerger.getMessages());

        TopoDS_Shape duplicateProcessedShape = duplicateVertexMerger.getResultShape();
        if (!duplicateProcessedShape.IsNull() &&
            duplicateProcessedShape.ShapeType() == originalShape.ShapeType() &&
            duplicateVertexMerger.isModified()) {
            workingShape = duplicateProcessedShape;
            duplicateMergedCount += containerDuplicateMergedCount;
            containerModified = true;
        }
        else if (duplicateVertexMerger.isModified()) {
            containerWarning = true;
            logger.container(Logger::Level::Warn, i, "process", logFields({
                Logger::textField("stage", "validate"),
                Logger::textField("action", "stage_rejected"),
                Logger::textField("algorithm", "MergeDuplicateVertex"),
                Logger::textField("fallback", "use_original_shape"),
                Logger::textField("reason", "result is null or changed shape type")
            }));
        }

        HealNarrowGap narrowGapHealer(narrowGapParams);
        narrowGapHealer.setInputShape(workingShape);
        const int containerNarrowGapRepairedCount = narrowGapHealer.perform();
        logger.containerMessages(i, "HealNarrowGap", narrowGapHealer.getMessages());

        TopoDS_Shape narrowGapProcessedShape = narrowGapHealer.getResultShape();
        if (!narrowGapProcessedShape.IsNull() &&
            narrowGapProcessedShape.ShapeType() == originalShape.ShapeType() &&
            narrowGapHealer.isModified()) {
            workingShape = narrowGapProcessedShape;
            narrowGapRepairedCount += containerNarrowGapRepairedCount;
            containerModified = true;
        }
        else if (narrowGapHealer.isModified()) {
            std::cerr << "Container " << i
                      << " narrow-gap result is invalid or changed shape type, keep previous shape." << std::endl;
            containerWarning = true;
            logger.container(Logger::Level::Warn, i, "process", logFields({
                Logger::textField("stage", "validate"),
                Logger::textField("action", "stage_rejected"),
                Logger::textField("algorithm", "HealNarrowGap"),
                Logger::textField("fallback", "keep_previous_shape"),
                Logger::textField("reason", "result is null or changed shape type")
            }));
        }

        if (!containerModified) {
            if (containerWarning) {
                ++warningCount;
            }
            logger.container(
                containerWarning ? Logger::Level::Warn : Logger::Level::Info,
                i,
                "end",
                logFields({
                    Logger::textField("status", containerWarning ? "failed" : "no_change"),
                    Logger::boolField("modified", false),
                    Logger::integerField("repair_count", 0),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, workingShape);
        ++replacedCount;
        if (containerWarning) {
            ++warningCount;
        }
        logger.container(
            containerWarning ? Logger::Level::Warn : Logger::Level::Info,
            i,
            "end",
            logFields({
                Logger::textField("status", containerWarning ? "partial_success" : "success"),
                Logger::boolField("modified", true),
                Logger::integerField(
                    "repair_count",
                    containerDuplicateMergedCount + containerNarrowGapRepairedCount),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after duplicate-vertex/narrow-gap repair: "
              << replacedCount << std::endl;
    std::cout << "Duplicate vertices merged: " << duplicateMergedCount
              << ", narrow gaps repaired: " << narrowGapRepairedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(warningCount > 0 ? Logger::Level::Warn : Logger::Level::Info, "end", logFields({
        Logger::textField("status", warningCount > 0 ? "partial_success" : "success"),
        Logger::integerField("container_count", static_cast<long long>(targets.size())),
        Logger::integerField("modified_containers", replacedCount),
        Logger::integerField("skipped_containers", skippedCount),
        Logger::integerField("warning_containers", warningCount),
        Logger::integerField("repair_count", duplicateMergedCount + narrowGapRepairedCount),
        Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
    }));

    return duplicateMergedCount + narrowGapRepairedCount;
}

int ApplicationInit::healNarrowGapsOnLowestContainersAndExport(const std::string& filePath,
                                                               const HealNarrowGap::Params& params)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path healLogPath = std::filesystem::path("assets") / "fixedNarrowGap" / "HealGapsLog.txt";
    Logger logger(healLogPath, "HealNarrowGap");
    reportLogOpenFailure(logger);
    logger.run(Logger::Level::Info, "start", logFields({
        Logger::pathField("output_step", filePath),
        Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
    }));

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot heal narrow gaps because the XCAF document is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int repairedCount = 0;
    int replacedCount = 0;
    int skippedCount = 0;
    int failedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "owner label is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "input shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        logger.container(Logger::Level::Info, i, "start", shapeTypeField(originalShape));

        HealNarrowGap narrowGapHealer(params);
        narrowGapHealer.setInputShape(originalShape);
        const int containerRepairedCount = narrowGapHealer.perform();
        logger.containerMessages(i, "HealNarrowGap", narrowGapHealer.getMessages());

        TopoDS_Shape processedShape = narrowGapHealer.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (!narrowGapHealer.isModified()) {
            logger.container(Logger::Level::Info, i, "end", logFields({
                Logger::textField("status", "no_change"), Logger::boolField("modified", false),
                Logger::integerField("repair_count", 0),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        repairedCount += containerRepairedCount;
        ++replacedCount;
        logger.container(Logger::Level::Info, i, "end", logFields({
            Logger::textField("status", "success"), Logger::boolField("modified", true),
            Logger::integerField("repair_count", containerRepairedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));

    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after narrow-gap heal: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(failedCount > 0 ? Logger::Level::Warn : Logger::Level::Info, "end", logFields({
        Logger::textField("status", failedCount > 0 ? "partial_success" : "success"),
        Logger::integerField("container_count", static_cast<long long>(targets.size())),
        Logger::integerField("modified_containers", replacedCount),
        Logger::integerField("skipped_containers", skippedCount),
        Logger::integerField("failed_containers", failedCount),
        Logger::integerField("repair_count", repairedCount),
        Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
    }));

    return repairedCount;
}

int ApplicationInit::healClosedGapsOnLowestContainersAndExport(
    const std::string& filePath,
    const HealClosedNarrowGap::Params& params)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path healLogPath =
        std::filesystem::path("assets") / "fixedClosedGap" / "HealClosedGapsLog.txt";
    Logger logger(healLogPath, "HealClosedNarrowGap");
    reportLogOpenFailure(logger);
    logger.run(Logger::Level::Info, "start", logFields({
        Logger::pathField("output_step", filePath),
        Logger::textField("log_level", Logger::levelName(logger.minimumLevel()))
    }));

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot heal closed gaps because the XCAF document is null." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int repairedCount = 0;
    int replacedCount = 0;
    int skippedCount = 0;
    int failedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "owner label is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            ++skippedCount;
            logger.container(Logger::Level::Warn, i, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "input shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        logger.container(Logger::Level::Info, i, "start", shapeTypeField(originalShape));

        HealClosedNarrowGap closedGapHealer(params);
        closedGapHealer.setInputShape(originalShape);
        const int containerRepairedCount = closedGapHealer.perform();
        logger.containerMessages(i, "HealClosedNarrowGap", closedGapHealer.getMessages());

        TopoDS_Shape processedShape = closedGapHealer.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            ++failedCount;
            logger.container(Logger::Level::Error, i, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (!closedGapHealer.isModified()) {
            logger.container(Logger::Level::Info, i, "end", logFields({
                Logger::textField("status", "no_change"), Logger::boolField("modified", false),
                Logger::integerField("repair_count", 0),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        repairedCount += containerRepairedCount;
        ++replacedCount;
        logger.container(Logger::Level::Info, i, "end", logFields({
            Logger::textField("status", "success"), Logger::boolField("modified", true),
            Logger::integerField("repair_count", containerRepairedCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after closed-gap heal: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(failedCount > 0 ? Logger::Level::Warn : Logger::Level::Info, "end", logFields({
        Logger::textField("status", failedCount > 0 ? "partial_success" : "success"),
        Logger::integerField("container_count", static_cast<long long>(targets.size())),
        Logger::integerField("modified_containers", replacedCount),
        Logger::integerField("skipped_containers", skippedCount),
        Logger::integerField("failed_containers", failedCount),
        Logger::integerField("repair_count", repairedCount),
        Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
    }));

    return repairedCount;
}

int ApplicationInit::runRepairPipelineOnLowestContainersAndExport(
    const std::string& filePath,
    const ShapeRepairPipeline& pipeline)
{
    const RepairClock::time_point runStart = RepairClock::now();
    const std::filesystem::path outputPath(filePath);
    const std::filesystem::path outputDirectory = outputPath.parent_path().empty()
        ? std::filesystem::path(".")
        : outputPath.parent_path();
    const std::filesystem::path pipelineLogPath =
        outputDirectory / "ShapeRepairPipelineLog.txt";
    Logger logger(pipelineLogPath, "ShapeRepairPipeline");
    reportLogOpenFailure(logger);
    logger.run(Logger::Level::Info, "start", logFields({
        Logger::pathField("output_step", filePath),
        Logger::textField("log_level", Logger::levelName(logger.minimumLevel())),
        Logger::integerField("configured_stage_count", static_cast<long long>(pipeline.getStages().size()))
    }));

    const std::vector<ShapeRepairPipeline::Stage>& configuredStages = pipeline.getStages();
    for (std::size_t stageIndex = 0; stageIndex < configuredStages.size(); ++stageIndex) {
        logger.run(Logger::Level::Info, "process", logFields({
            Logger::textField("stage", "config"),
            Logger::textField("action", "stage_configured"),
            Logger::integerField("stage_index", static_cast<long long>(stageIndex)),
            Logger::textField(
                "algorithm",
                configuredStages[stageIndex].name.empty() ? "UnnamedStage" : configuredStages[stageIndex].name),
            Logger::boolField("enabled", configuredStages[stageIndex].enabled)
        }));
    }

    if (pipeline.empty()) {
        std::cerr << "Error: Cannot run an empty repair pipeline." << std::endl;
        logRunFailure(logger, runStart, "config", "repair pipeline is empty");
        return 0;
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot run repair pipeline because the XCAF document is null."
                  << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF document is null");
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        logRunFailure(logger, runStart, "input_validate", "XCAF shape tool is null");
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        logRunFailure(logger, runStart, "container_collect", "no editable solid/shell container");
        return 0;
    }

    int totalRepairCount = 0;
    int replacedContainerCount = 0;
    int skippedContainerCount = 0;
    int failedContainerCount = 0;
    for (std::size_t containerIndex = 0; containerIndex < targets.size(); ++containerIndex) {
        const RepairClock::time_point containerStart = RepairClock::now();
        const EditableContainerTarget& target = targets[containerIndex];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << containerIndex << " has no owner label, skip."
                      << std::endl;
            ++skippedContainerCount;
            logger.container(Logger::Level::Warn, containerIndex, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "owner label is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull()
            ? target.originalShape
            : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << containerIndex << " has null shape, skip."
                      << std::endl;
            ++skippedContainerCount;
            logger.container(Logger::Level::Warn, containerIndex, "end", logFields({
                Logger::textField("status", "skipped"),
                Logger::textField("reason", "input shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        logger.container(Logger::Level::Info, containerIndex, "start", shapeTypeField(originalShape));

        const ShapeRepairPipeline::Result pipelineResult = pipeline.perform(originalShape);

        // 此时管线已经执行完所有启用阶段。按执行顺序统一写出阶段报告，
        // 让各算法无需感知日志路径、Container 编号和当前日志等级。
        for (std::size_t stageIndex = 0;
             stageIndex < pipelineResult.stageReports.size();
             ++stageIndex) {
            const ShapeRepairPipeline::StageReport& report = pipelineResult.stageReports[stageIndex];
            const char* stageStatus = "skipped";
            Logger::Level stageLevel = Logger::Level::Info;
            if (!report.enabled) {
                stageStatus = "disabled";
            }
            else if (!report.rejectionReason.empty()) {
                stageStatus = "rejected";
                stageLevel = Logger::Level::Warn;
            }
            else if (report.accepted && report.modified) {
                stageStatus = "accepted";
            }
            else if (report.accepted) {
                stageStatus = "no_change";
            }

            logger.container(stageLevel, containerIndex, "process", logFields({
                Logger::textField("stage", "repair"),
                Logger::textField("action", "stage_completed"),
                Logger::integerField("stage_index", static_cast<long long>(stageIndex)),
                Logger::textField("algorithm", report.name),
                Logger::textField("status", stageStatus),
                Logger::boolField("executed", report.executed),
                Logger::boolField("modified", report.modified),
                Logger::integerField("repair_count", report.repairCount),
                Logger::textField("input_shape_type", TopologyUtils::ShapeTypeName(report.inputShapeType)),
                Logger::textField("output_shape_type", TopologyUtils::ShapeTypeName(report.outputShapeType))
            }));

            if (!report.rejectionReason.empty()) {
                // 阶段结果被拒绝并不代表整个管线立即失败：该结果不会提交，
                // 后续阶段仍可基于上一份已通过校验的 Shape 继续执行。
                logger.containerText(Logger::Level::Warn, containerIndex, "process", logFields({
                    Logger::textField("stage", "validate"),
                    Logger::textField("action", "stage_rejected"),
                    Logger::textField("algorithm", report.name)
                }), "reason", report.rejectionReason);
            }
            logger.containerMessages(containerIndex, report.name, report.messages);
        }

        logger.containerMessages(containerIndex, "ShapeRepairPipeline", pipelineResult.messages);

        if (!pipelineResult.modified) {
            // “正常但没有修改”和“最终校验失败后回滚”都会保留原 Container，
            // 但后者属于操作失败，必须用 ERROR 和 rolled_back=true 明确区分。
            if (pipelineResult.rolledBack) {
                ++failedContainerCount;
                logger.container(Logger::Level::Error, containerIndex, "end", logFields({
                    Logger::textField("status", "failed"),
                    Logger::boolField("modified", false),
                    Logger::boolField("rolled_back", true),
                    Logger::textField("action", "keep_original"),
                    Logger::textField("reason", "final pipeline validation rejected the result"),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            }
            else {
                logger.container(Logger::Level::Info, containerIndex, "end", logFields({
                    Logger::textField("status", "no_change"),
                    Logger::boolField("modified", false),
                    Logger::integerField("repair_count", 0),
                    Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
                }));
            }
            continue;
        }

        if (pipelineResult.resultShape.IsNull()) {
            std::cerr << "Container " << containerIndex
                      << " pipeline produced a null result, keep original shape."
                      << std::endl;
            ++failedContainerCount;
            logger.container(Logger::Level::Error, containerIndex, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape is null"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        if (pipelineResult.resultShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << containerIndex
                      << " pipeline changed shape type, keep original shape."
                      << std::endl;
            ++failedContainerCount;
            logger.container(Logger::Level::Error, containerIndex, "end", logFields({
                Logger::textField("status", "failed"), Logger::boolField("modified", false),
                Logger::textField("action", "keep_original"),
                Logger::textField("reason", "result shape type changed"),
                Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
            }));
            continue;
        }

        // 整条链完成后只回写一次，避免中间结果使 XCAF reference 反复失效。
        shapeTool->SetShape(target.ownerLabel, pipelineResult.resultShape);
        totalRepairCount += pipelineResult.totalRepairCount;
        ++replacedContainerCount;
        logger.container(Logger::Level::Info, containerIndex, "end", logFields({
            Logger::textField("status", "success"),
            Logger::boolField("modified", true),
            Logger::boolField("rolled_back", false),
            Logger::integerField("repair_count", pipelineResult.totalRepairCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(containerStart))
        }));
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after repair pipeline: "
              << replacedContainerCount << std::endl;
    std::cout << "Repair pipeline committed operation count: "
              << totalRepairCount << std::endl;

    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        logRunFailure(logger, runStart, "export", "failed to write STEP document");
        return -1;
    }

    logger.run(
        failedContainerCount > 0 ? Logger::Level::Warn : Logger::Level::Info,
        "end",
        logFields({
            Logger::textField("status", failedContainerCount > 0 ? "partial_success" : "success"),
            Logger::integerField("container_count", static_cast<long long>(targets.size())),
            Logger::integerField("modified_containers", replacedContainerCount),
            Logger::integerField("skipped_containers", skippedContainerCount),
            Logger::integerField("failed_containers", failedContainerCount),
            Logger::integerField("repair_count", totalRepairCount),
            Logger::integerField("duration_ms", elapsedMilliseconds(runStart))
        }));

    return totalRepairCount;
}

bool ApplicationInit::writeXcafDocumentStep(const Handle(TDocStd_Document)& doc,
                                            const std::string& filePath,
                                            bool preserveAssembly) const
{
    if (doc.IsNull()) {
        std::cerr << "Error: The XCAF document to export is null." << std::endl;
        return false;
    }

    std::cout << "Writing XCAF STEP file to: " << filePath << "......" << std::endl;
    Interface_Static::SetIVal("write.step.assembly", preserveAssembly ? 1 : 0);
    setOutputStepUnitFromInput();

    STEPCAFControl_Writer writer;

    ScopedOccInfoSilencer occInfoSilencer;
    if (!writer.Transfer(doc, STEPControl_AsIs)) {
        std::cerr << "Error: Failed to transfer XCAF document to STEP writer." << std::endl;
        return false;
    }

    if (writer.Write(filePath.c_str()) == IFSelect_RetDone) {
        std::cout << "XCAF STEP file written successfully." << std::endl;
        return true;
    }

    std::cerr << "Error: Failed to write XCAF STEP file." << std::endl;
    return false;
}

bool ApplicationInit::transferAndCacheTopLevelShapes(Standard_Boolean isTransferred, const char* fileType)
{
    if (!isTransferred) {
        std::cerr << "Error: Failed to transfer " << fileType << " file to document." << std::endl;
        clearShapeCache();
        return false;
    }

    std::cout << fileType << " file transferred to document successfully." << std::endl;
    refreshShapeCache();
    return true;
}

void ApplicationInit::refreshShapeCache()
{
    clearShapeCache();

    if (aDoc_.IsNull()) {
        return;
    }
    
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        return;
    }
    TDF_LabelSequence freeShapes;
	TDF_LabelSequence topLevelShapes;
    TopTools_MapOfShape shapeMap;
    shapeTool->GetFreeShapes(freeShapes);
	std::cout << "Free shapes from document, count: " << freeShapes.Length() << std::endl;
    shapeTool->GetShapes(topLevelShapes);
    std::cout << "Top-level shapes from document, count: " << topLevelShapes.Length() << std::endl;

    for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
        const TopoDS_Shape aShape = shapeTool->GetShape(freeShapes.Value(i));
        if (aShape.IsNull() || !shapeMap.Add(aShape)) {
            continue;
        }
        aShapeVec_.push_back(aShape);
        if (aShape.ShapeType() == TopAbs_SOLID) {
            aSolidVec_.push_back(TopoDS::Solid(aShape));
        }
        else if (aShape.ShapeType() == TopAbs_SHELL) {
            aShellVec_.push_back(TopoDS::Shell(aShape));
        }
    }
}

void ApplicationInit::clearShapeCache()
{
    aShapeVec_.clear();
    aSolidVec_.clear();
    aShellVec_.clear();
}

bool ApplicationInit::shapeContainsSolidOrShell(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return false;
    }

    if (shape.ShapeType() == TopAbs_SOLID || shape.ShapeType() == TopAbs_SHELL) {
        return true;
    }

    TopExp_Explorer solidExplorer(shape, TopAbs_SOLID);
    if (solidExplorer.More()) {
        return true;
    }

    TopExp_Explorer shellExplorer(shape, TopAbs_SHELL);
    return shellExplorer.More();
}

void ApplicationInit::setOutputStepUnitFromInput() const
{
    Interface_Static::SetCVal("write.step.unit", inputStepUnit_.c_str());
    std::cout << "Output STEP length unit: " << inputStepUnit_ << std::endl;
}

std::vector<TopoDS_Shape> ApplicationInit::getAllShapes() const
{
    return aShapeVec_;
}

std::vector<TopoDS_Solid> ApplicationInit::getAllSolids() const
{
    return aSolidVec_;
}

std::vector<TopoDS_Shell> ApplicationInit::getAllShells() const
{
    return aShellVec_;
}

std::vector<TopoDS_Edge> ApplicationInit::getAllEdges() const
{
    std::vector<TopoDS_Edge> edgeVec;
    TopTools_MapOfShape edgeMap;

    for (const auto& shape : aShapeVec_) {
        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            const TopoDS_Shape edgeShape = exp.Current();
            if (!edgeMap.Add(edgeShape)) {
                continue;
            }
            edgeVec.push_back(TopoDS::Edge(edgeShape));
        }
    }

    return edgeVec;
}
