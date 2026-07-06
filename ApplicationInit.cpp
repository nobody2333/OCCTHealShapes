#include "ApplicationInit.h"

#include <cctype>
#include <filesystem>
#include <fstream>
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

#include "TopologyUtils.h"

namespace {
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
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexStepShape" / "MergeDuplicateVertexStepShapeLog.txt";
    std::filesystem::create_directories(mergeLogPath.parent_path());
    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "MergeDuplicateVertex STEP shape log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open MergeDuplicateVertex STEP shape log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: Cannot merge duplicate vertices because imported STEP shape is null." << std::endl;
        return 0;
    }

    MergeDuplicateVertex duplicateVertexMerger(params);
    duplicateVertexMerger.setInputShape(importedStepShape_);
    const int mergedCount = duplicateVertexMerger.perform();

    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::app);
        if (logFile) {
            logFile << "[MergeDuplicateVertex][StepShape]"
                    << " mergedCount=" << mergedCount
                    << ", modified=" << (duplicateVertexMerger.isModified() ? "true" : "false")
                    << ", inputShapeType=" << static_cast<int>(importedStepShape_.ShapeType())
                    << std::endl;

            const std::vector<std::string>& messages = duplicateVertexMerger.getMessages();
            if (messages.empty()) {
                logFile << "  (no messages)" << std::endl;
            }
            for (const std::string& message : messages) {
                logFile << "  - " << message << std::endl;
            }
            logFile << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot append MergeDuplicateVertex STEP shape log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    TopoDS_Shape processedShape = duplicateVertexMerger.getResultShape();
    if (processedShape.IsNull()) {
        std::cerr << "Error: MergeDuplicateVertex produced a null STEP shape result." << std::endl;
        return 0;
    }

    if (processedShape.ShapeType() != importedStepShape_.ShapeType()) {
        std::cerr << "Error: MergeDuplicateVertex changed STEP shape type from "
                  << static_cast<int>(importedStepShape_.ShapeType()) << " to "
                  << static_cast<int>(processedShape.ShapeType()) << "." << std::endl;
        return 0;
    }

    importedStepShape_ = processedShape;

    if (!exportStepShape(importedStepShape_, filePath)) {
        return -1;
    }

    return mergedCount;
}

int ApplicationInit::sewShellsOnStepShapeAndExport(
    const std::string& filePath,
    const SewShells::Params& params)
{
    const std::filesystem::path sewLogPath =
        std::filesystem::path("assets") / "fixedSewShellsStepShape" / "SewShellsStepShapeLog.txt";
    std::filesystem::create_directories(sewLogPath.parent_path());
    {
        std::ofstream logFile(sewLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "SewShells STEP shape log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open SewShells STEP shape log file: "
                      << sewLogPath.string() << std::endl;
        }
    }

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: Cannot sew shells because imported STEP shape is null." << std::endl;
        return 0;
    }

    SewShells shellSewer(params);
    shellSewer.setInputShape(importedStepShape_);
    const int sewedGroupCount = shellSewer.perform();

    {
        std::ofstream logFile(sewLogPath, std::ios::out | std::ios::app);
        if (logFile) {
            logFile << "[SewShells][StepShape]"
                    << " sewedGroupCount=" << sewedGroupCount
                    << ", modified=" << (shellSewer.isModified() ? "true" : "false")
                    << ", inputShapeType=" << TopologyUtils::ShapeTypeName(importedStepShape_)
                    << std::endl;

            const std::vector<std::string>& messages = shellSewer.getMessages();
            if (messages.empty()) {
                logFile << "  (no messages)" << std::endl;
            }
            for (const std::string& message : messages) {
                logFile << "  - " << message << std::endl;
            }
            logFile << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot append SewShells STEP shape log file: "
                      << sewLogPath.string() << std::endl;
        }
    }

    TopoDS_Shape processedShape = shellSewer.getResultShape();
    if (processedShape.IsNull()) {
        std::cerr << "Error: SewShells produced a null STEP shape result." << std::endl;
        return 0;
    }

    if (processedShape.ShapeType() != importedStepShape_.ShapeType()) {
        std::cerr << "Error: SewShells changed STEP shape type from "
                  << TopologyUtils::ShapeTypeName(importedStepShape_) << " to "
                  << TopologyUtils::ShapeTypeName(processedShape) << "." << std::endl;
        return 0;
    }

    importedStepShape_ = processedShape;

    if (!exportStepShape(importedStepShape_, filePath)) {
        return -1;
    }

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
    const std::filesystem::path mergeLogPath = std::filesystem::path("assets") / "fixedShortEdge" / "MergeShortEdgeLog.txt";
    std::filesystem::create_directories(mergeLogPath.parent_path());
    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "MergeShortEdge log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open MergeShortEdge log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot merge short edges because the XCAF document is null." << std::endl;
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    //collectEditableContainerTargets(targets);
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        return 0;
    }

    int mergedCount = 0;
    int replacedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            continue;
        }

        // 重新读取 label 当前 shape，避免前面容器写回后导致装配派生的
        // compound 已刷新，而这里仍使用旧缓存。
        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            continue;
        }

        MergeShortEdge shortEdgeMerger;
        if (useManualLengthThreshold) {
            shortEdgeMerger.setLengthThreshold(manualLengthThreshold);
        }
        shortEdgeMerger.setAngleThresholdDeg(angleThresholdDeg);
        shortEdgeMerger.setInputShape(originalShape);
        const int containerMergedCount = shortEdgeMerger.performMerge();

        {
            std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[MergeShortEdge] Container " << i
                        << ", mergedCount=" << containerMergedCount
                        << ", inputShapeType=" << static_cast<int>(originalShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = shortEdgeMerger.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
                logFile << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append MergeShortEdge log file: "
                          << mergeLogPath.string() << std::endl;
            }
        }

        TopoDS_Shape processedShape = shortEdgeMerger.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            continue;
        }

        // 整个叶子容器只写回一次。这样会牺牲该局部装配 label 内部的
        // 细粒度装配语义，但能保留上层 STEP/XCAF 父级关系，避免 Cadfix
        // 将 shell 下的 face 识别成 parentless face。
        shapeTool->SetShape(target.ownerLabel, processedShape);
        mergedCount += containerMergedCount;
        ++replacedCount;
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after short-edge merge: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return mergedCount;
}

int ApplicationInit::mergeDuplicateVerticesOnLowestContainersAndExport(
    const std::string& filePath,
    const MergeDuplicateVertex::Params& params)
{
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertex" / "MergeDuplicateVertexLog.txt";
    std::filesystem::create_directories(mergeLogPath.parent_path());
    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "MergeDuplicateVertex log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open MergeDuplicateVertex log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot merge duplicate vertices because the XCAF document is null." << std::endl;
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        return 0;
    }

    int mergedCount = 0;
    int replacedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            continue;
        }

        MergeDuplicateVertex duplicateVertexMerger(params);
        duplicateVertexMerger.setInputShape(originalShape);
        const int containerMergedCount = duplicateVertexMerger.perform();

        {
            std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[MergeDuplicateVertex] Container " << i
                        << ", mergedCount=" << containerMergedCount
                        << ", modified=" << (duplicateVertexMerger.isModified() ? "true" : "false")
                        << ", inputShapeType=" << static_cast<int>(originalShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = duplicateVertexMerger.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
                logFile << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append MergeDuplicateVertex log file: "
                          << mergeLogPath.string() << std::endl;
            }
        }

        TopoDS_Shape processedShape = duplicateVertexMerger.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            continue;
        }

        if (!duplicateVertexMerger.isModified()) {
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        mergedCount += containerMergedCount;
        ++replacedCount;
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after duplicate-vertex merge: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return mergedCount;
}

int ApplicationInit::mergeDuplicateVerticesAndHealNarrowGapsOnLowestContainersAndExport(
    const std::string& filePath,
    const MergeDuplicateVertex::Params& duplicateVertexParams,
    const HealNarrowGap::Params& narrowGapParams)
{
    const std::filesystem::path chainLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexNarrowGap" / "MergeDuplicateVertexHealNarrowGapLog.txt";
    std::filesystem::create_directories(chainLogPath.parent_path());
    {
        std::ofstream logFile(chainLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "MergeDuplicateVertex + HealNarrowGap log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open duplicate-vertex/narrow-gap log file: "
                      << chainLogPath.string() << std::endl;
        }
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot run duplicate-vertex/narrow-gap repair because the XCAF document is null." << std::endl;
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        return 0;
    }

    int duplicateMergedCount = 0;
    int narrowGapRepairedCount = 0;
    int replacedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            continue;
        }

        TopoDS_Shape workingShape = originalShape;
        bool containerModified = false;

        MergeDuplicateVertex duplicateVertexMerger(duplicateVertexParams);
        duplicateVertexMerger.setInputShape(workingShape);
        const int containerDuplicateMergedCount = duplicateVertexMerger.perform();

        {
            std::ofstream logFile(chainLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[MergeDuplicateVertex] Container " << i
                        << ", mergedCount=" << containerDuplicateMergedCount
                        << ", modified=" << (duplicateVertexMerger.isModified() ? "true" : "false")
                        << ", inputShapeType=" << static_cast<int>(workingShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = duplicateVertexMerger.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
            }
        }

        TopoDS_Shape duplicateProcessedShape = duplicateVertexMerger.getResultShape();
        if (!duplicateProcessedShape.IsNull() &&
            duplicateProcessedShape.ShapeType() == originalShape.ShapeType() &&
            duplicateVertexMerger.isModified()) {
            workingShape = duplicateProcessedShape;
            duplicateMergedCount += containerDuplicateMergedCount;
            containerModified = true;
        }
        else if (duplicateVertexMerger.isModified()) {
            std::ofstream logFile(chainLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "  - Warning: duplicate-vertex result is invalid or changed shape type; "
                        << "HealNarrowGap will use original container shape." << std::endl;
            }
        }

        HealNarrowGap narrowGapHealer(narrowGapParams);
        narrowGapHealer.setInputShape(workingShape);
        const int containerNarrowGapRepairedCount = narrowGapHealer.perform();

        {
            std::ofstream logFile(chainLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[HealNarrowGap] Container " << i
                        << ", repairedCount=" << containerNarrowGapRepairedCount
                        << ", modified=" << (narrowGapHealer.isModified() ? "true" : "false")
                        << ", inputShapeType=" << static_cast<int>(workingShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = narrowGapHealer.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
                logFile << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append duplicate-vertex/narrow-gap log file: "
                          << chainLogPath.string() << std::endl;
            }
        }

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
        }

        if (!containerModified) {
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, workingShape);
        ++replacedCount;
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after duplicate-vertex/narrow-gap repair: "
              << replacedCount << std::endl;
    std::cout << "Duplicate vertices merged: " << duplicateMergedCount
              << ", narrow gaps repaired: " << narrowGapRepairedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return duplicateMergedCount + narrowGapRepairedCount;
}

int ApplicationInit::healNarrowGapsOnLowestContainersAndExport(const std::string& filePath,
                                                               const HealNarrowGap::Params& params)
{
    const std::filesystem::path healLogPath = std::filesystem::path("assets") / "fixedNarrowGap" / "HealGapsLog.txt";
    std::filesystem::create_directories(healLogPath.parent_path());
    {
        std::ofstream logFile(healLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "HealNarrowGap log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open HealNarrowGap log file: "
                      << healLogPath.string() << std::endl;
        }
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot heal narrow gaps because the XCAF document is null." << std::endl;
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        return 0;
    }

    int repairedCount = 0;
    int replacedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            continue;
        }

        HealNarrowGap narrowGapHealer(params);
        narrowGapHealer.setInputShape(originalShape);
        const int containerRepairedCount = narrowGapHealer.perform();

        {
            std::ofstream logFile(healLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[HealNarrowGap] Container " << i
                        << ", repairedCount=" << containerRepairedCount
                        << ", modified=" << (narrowGapHealer.isModified() ? "true" : "false")
					    << ", inputShapeType=" << static_cast<int>(originalShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = narrowGapHealer.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
                logFile << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append HealNarrowGap log file: "
                          << healLogPath.string() << std::endl;
            }
        }

        TopoDS_Shape processedShape = narrowGapHealer.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            continue;
        }

        if (!narrowGapHealer.isModified()) {
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        repairedCount += containerRepairedCount;
        ++replacedCount;

    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after narrow-gap heal: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return repairedCount;
}

int ApplicationInit::healClosedGapsOnLowestContainersAndExport(
    const std::string& filePath,
    const HealClosedNarrowGap::Params& params)
{
    const std::filesystem::path healLogPath =
        std::filesystem::path("assets") / "fixedClosedGap" / "HealClosedGapsLog.txt";
    std::filesystem::create_directories(healLogPath.parent_path());
    {
        std::ofstream logFile(healLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "HealClosedNarrowGap log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open HealClosedNarrowGap log file: "
                      << healLogPath.string() << std::endl;
        }
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot heal closed gaps because the XCAF document is null." << std::endl;
        return 0;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        std::cerr << "Error: Failed to get XCAF shape tool." << std::endl;
        return 0;
    }

    std::vector<EditableContainerTarget> targets;
    collectTopLevelTargets(targets);
    if (targets.empty()) {
        std::cerr << "Error: Cannot find any editable solid/shell container." << std::endl;
        return 0;
    }

    int repairedCount = 0;
    int replacedCount = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const EditableContainerTarget& target = targets[i];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << i << " has no owner label, skip." << std::endl;
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull() ? target.originalShape : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << i << " has null shape, skip." << std::endl;
            continue;
        }

        HealClosedNarrowGap closedGapHealer(params);
        closedGapHealer.setInputShape(originalShape);
        const int containerRepairedCount = closedGapHealer.perform();

        {
            std::ofstream logFile(healLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[HealClosedNarrowGap] Container " << i
                        << ", repairedCount=" << containerRepairedCount
                        << ", modified=" << (closedGapHealer.isModified() ? "true" : "false")
                        << ", inputShapeType=" << static_cast<int>(originalShape.ShapeType())
                        << std::endl;

                const std::vector<std::string>& messages = closedGapHealer.getMessages();
                if (messages.empty()) {
                    logFile << "  (no messages)" << std::endl;
                }
                for (const std::string& message : messages) {
                    logFile << "  - " << message << std::endl;
                }
                logFile << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append HealClosedNarrowGap log file: "
                          << healLogPath.string() << std::endl;
            }
        }

        TopoDS_Shape processedShape = closedGapHealer.getResultShape();
        if (processedShape.IsNull()) {
            std::cerr << "Container " << i << " produced a null result, keep original shape." << std::endl;
            continue;
        }

        if (processedShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << i << " changed shape type from "
                      << static_cast<int>(originalShape.ShapeType()) << " to "
                      << static_cast<int>(processedShape.ShapeType())
                      << ", keep original shape." << std::endl;
            continue;
        }

        if (!closedGapHealer.isModified()) {
            continue;
        }

        shapeTool->SetShape(target.ownerLabel, processedShape);
        repairedCount += containerRepairedCount;
        ++replacedCount;
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after closed-gap heal: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return repairedCount;
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
