#include "ApplicationInit.h"

#include <cctype>
#include <algorithm>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include "ShapeRepairStages.h"
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
#include "RepairDiagnostics.h"

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
    const ShapePreProcess::Params& params)
{
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexStepShape" / "MergeDuplicateVertexStepShapeLog.txt";
    std::filesystem::create_directories(mergeLogPath.parent_path());
    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "Shape preprocessing (STEP shape)" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open preprocessing STEP log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    if (importedStepShape_.IsNull()) {
        std::cerr << "Error: Cannot merge duplicate vertices because imported STEP shape is null." << std::endl;
        return 0;
    }

    ShapePreProcess duplicateVertexMerger(params);
    duplicateVertexMerger.setInputShape(importedStepShape_);
    const int mergedCount = duplicateVertexMerger.perform();

    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::app);
        if (logFile) {
            logFile << "[ShapePreProcess][StepShape]"
                    << " mergedPairs=" << duplicateVertexMerger.getMergedCount()
                    << ", shortEdgeCases=" << duplicateVertexMerger.getCollapsedBoundaryCount()
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
            std::cerr << "Warning: Cannot append preprocessing STEP log file: "
                      << mergeLogPath.string() << std::endl;
        }
    }

    TopoDS_Shape processedShape = duplicateVertexMerger.getResultShape();
    if (processedShape.IsNull()) {
        std::cerr << "Error: ShapePreProcess produced a null STEP shape result." << std::endl;
        return 0;
    }

    if (processedShape.ShapeType() != importedStepShape_.ShapeType()) {
        std::cerr << "Error: ShapePreProcess changed STEP shape type from "
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
    const ShapePreProcess::Params& params)
{
    const std::filesystem::path mergeLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertex" / "MergeDuplicateVertexLog.txt";
    std::filesystem::create_directories(mergeLogPath.parent_path());
    {
        std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "Shape preprocessing" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open preprocessing log file: "
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

        ShapePreProcess duplicateVertexMerger(params);
        duplicateVertexMerger.setInputShape(originalShape);
        const int containerMergedCount = duplicateVertexMerger.perform();

        {
            std::ofstream logFile(mergeLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[ShapePreProcess] Container " << i
                        << ", mergedPairs=" << duplicateVertexMerger.getMergedCount()
                        << ", shortEdgeCases=" << duplicateVertexMerger.getCollapsedBoundaryCount()
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
                std::cerr << "Warning: Cannot append preprocessing log file: "
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
    std::cout << "Replaced lowest containers after preprocessing: " << replacedCount << std::endl;
    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return mergedCount;
}

int ApplicationInit::mergeDuplicateVerticesAndHealNarrowGapsOnLowestContainersAndExport(
    const std::string& filePath,
    const ShapePreProcess::Params& duplicateVertexParams,
    const HealNarrowGap::Params& narrowGapParams)
{
    const std::filesystem::path chainLogPath =
        std::filesystem::path("assets") / "fixedDuplicateVertexNarrowGap" / "MergeDuplicateVertexHealNarrowGapLog.txt";
    std::filesystem::create_directories(chainLogPath.parent_path());
    {
        std::ofstream logFile(chainLogPath, std::ios::out | std::ios::trunc);
        if (logFile) {
            logFile << "Shape preprocessing + HealNarrowGap" << std::endl;
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

        ShapePreProcess duplicateVertexMerger(duplicateVertexParams);
        duplicateVertexMerger.setInputShape(workingShape);
        const int containerDuplicateMergedCount = duplicateVertexMerger.perform();

        {
            std::ofstream logFile(chainLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[ShapePreProcess] Container " << i
                        << ", mergedPairs=" << duplicateVertexMerger.getMergedCount()
                        << ", shortEdgeCases=" << duplicateVertexMerger.getCollapsedBoundaryCount()
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
                logFile << "  - Warning: preprocessing result is invalid or changed shape type; "
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
                    if (message.empty()) {
                        logFile << std::endl;
                        continue;
                    }
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
    std::cout << "Replaced lowest containers after preprocessing/narrow-gap repair: "
              << replacedCount << std::endl;
    std::cout << "Preprocessing operations: " << duplicateMergedCount
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
                    if (message.empty()) {
                        logFile << std::endl;
                        continue;
                    }
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
                    if (message.empty()) {
                        logFile << std::endl;
                        continue;
                    }
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

int ApplicationInit::runRepairPipelineOnLowestContainersAndExport(
    const std::string& filePath,
    const ShapeRepairPipeline& pipeline)
{
    const std::filesystem::path outputPath(filePath);
    const std::filesystem::path outputDirectory = outputPath.parent_path().empty()
        ? std::filesystem::path(".") : outputPath.parent_path();
    const std::filesystem::path pipelineLogPath =
        outputDirectory / "ShapeRepairPipelineLog.txt";
    std::filesystem::create_directories(outputDirectory);

    {
        std::ofstream logFile(pipelineLogPath, std::ios::trunc);
        if (logFile) {
            logFile << "ShapeRepairPipeline log" << std::endl;
            logFile << "Output STEP: " << filePath << std::endl;
            logFile << "Configured stages:" << std::endl;
            const std::vector<ShapeRepairPipeline::Stage>& stages = pipeline.getStages();
            if (stages.empty()) {
                logFile << "  (none)" << std::endl;
            }
            for (std::size_t stageIndex = 0; stageIndex < stages.size(); ++stageIndex) {
                logFile << "  " << stageIndex << ": "
                        << (stages[stageIndex].name.empty() ? "UnnamedStage" : stages[stageIndex].name)
                        << ", enabled=" << (stages[stageIndex].enabled ? "true" : "false")
                        << std::endl;
            }
            logFile << std::endl;
        }
        else {
            std::cerr << "Warning: Cannot open repair pipeline log file: "
                      << pipelineLogPath.string() << std::endl;
        }
    }

    if (pipeline.empty()) {
        std::cerr << "Error: Cannot run an empty repair pipeline." << std::endl;
        return 0;
    }

    if (aDoc_.IsNull()) {
        std::cerr << "Error: Cannot run repair pipeline because the XCAF document is null."
                  << std::endl;
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

    int totalRepairCount = 0;
    int replacedContainerCount = 0;
    for (std::size_t containerIndex = 0; containerIndex < targets.size(); ++containerIndex) {
        const EditableContainerTarget& target = targets[containerIndex];
        if (target.ownerLabel.IsNull()) {
            std::cerr << "Container " << containerIndex << " has no owner label, skip."
                      << std::endl;
            continue;
        }

        const TopoDS_Shape currentShape = shapeTool->GetShape(target.ownerLabel);
        const TopoDS_Shape originalShape = currentShape.IsNull()
            ? target.originalShape
            : currentShape;
        if (originalShape.IsNull()) {
            std::cerr << "Container " << containerIndex << " has null shape, skip."
                      << std::endl;
            continue;
        }

        const ShapeRepairPipeline::Result pipelineResult = pipeline.perform(originalShape);

        {
            std::ofstream logFile(pipelineLogPath, std::ios::out | std::ios::app);
            if (logFile) {
                logFile << "[Container " << containerIndex << "]"
                        << ", inputShapeType=" << TopologyUtils::ShapeTypeName(originalShape)
                        << std::endl << std::endl;

                for (std::size_t stageIndex = 0;
                     stageIndex < pipelineResult.stageReports.size();
                     ++stageIndex) {
                    const ShapeRepairPipeline::StageReport& report =
                        pipelineResult.stageReports[stageIndex];

                    const char* status = "skipped";
                    if (!report.enabled) {
                        status = "disabled";
                    }
                    else if (!report.rejectionReason.empty()) {
                        status = "rejected";
                    }
                    else if (report.accepted && report.modified) {
                        status = "accepted";
                    }
                    else if (report.accepted) {
                        status = "no-change";
                    }

                    logFile << "  [Stage " << stageIndex << ": " << report.name << "]"
                            << std::endl;
                    logFile << "    status=" << status
                            << ", executed=" << (report.executed ? "true" : "false")
                            << ", modified=" << (report.modified ? "true" : "false")
                            << ", repairCount=" << report.repairCount
                            << ", inputShapeType=" << TopologyUtils::ShapeTypeName(report.inputShapeType)
                            << ", outputShapeType=" << TopologyUtils::ShapeTypeName(report.outputShapeType)
                            << std::endl;

                    if (!report.rejectionReason.empty()) {
                        logFile << "    rejectionReason:" << std::endl
                                << RepairDiagnostics::IndentMultiline(
                                       RepairDiagnostics::MakeFailureReasonReadable(
                                           report.rejectionReason),
                                       "      ")
                                << std::endl;
                    }

                    if (report.messages.empty()) {
                        logFile << "    (no messages)" << std::endl;
                    }
                    for (const std::string& message : report.messages) {
                        if (message.empty()) {
                            logFile << std::endl;
                            continue;
                        }
                        logFile << RepairDiagnostics::IndentMultiline(
                                       "- " + message,
                                       "    ")
                                << std::endl;
                    }
                    logFile << std::endl;
                }

                if (!pipelineResult.messages.empty()) {
                    logFile << "  [Pipeline messages]" << std::endl;
                    for (const std::string& message : pipelineResult.messages) {
                        logFile << RepairDiagnostics::IndentMultiline(
                                       "- " + message,
                                       "    ")
                                << std::endl;
                    }
                    logFile << std::endl;
                }

                logFile << "  [Pipeline result]"
                        << ", modified=" << (pipelineResult.modified ? "true" : "false")
                        << ", rolledBack=" << (pipelineResult.rolledBack ? "true" : "false")
                        << ", totalRepairCount=" << pipelineResult.totalRepairCount
                        << std::endl << std::endl;
            }
            else {
                std::cerr << "Warning: Cannot append repair pipeline log file: "
                          << pipelineLogPath.string() << std::endl;
            }
        }

        if (!pipelineResult.modified) {
            continue;
        }

        if (pipelineResult.resultShape.IsNull()) {
            std::cerr << "Container " << containerIndex
                      << " pipeline produced a null result, keep original shape."
                      << std::endl;
            continue;
        }

        if (pipelineResult.resultShape.ShapeType() != originalShape.ShapeType()) {
            std::cerr << "Container " << containerIndex
                      << " pipeline changed shape type, keep original shape."
                      << std::endl;
            continue;
        }

        // 整条链完成后只回写一次，避免中间结果使 XCAF reference 反复失效。
        shapeTool->SetShape(target.ownerLabel, pipelineResult.resultShape);
        totalRepairCount += pipelineResult.totalRepairCount;
        ++replacedContainerCount;
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after repair pipeline: "
              << replacedContainerCount << std::endl;
    std::cout << "Repair pipeline committed operation count: "
              << totalRepairCount << std::endl;

    if (!writeXcafDocumentStep(aDoc_, filePath, true)) {
        return -1;
    }

    return totalRepairCount;
}

int ApplicationInit::runGroupedRepairPipelineAndExport(
    const std::string& filePath, const ShapeRepairPipeline& pipeline)
{
    if (aDoc_.IsNull() || pipeline.empty()) return -1;
    const auto shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) return -1;
    const auto directory = std::filesystem::path(filePath).parent_path();
    if (!directory.empty()) std::filesystem::create_directories(directory);
    std::ofstream log(directory / "ShapeRepairPipelineLog.txt", std::ios::trunc);
    if (!log) return -1;

    // GetFreeShapes 只取实例根；GetShape 返回装配的实际位置化几何。
    // 不再同时遍历零件定义，避免把定义和实例重复输出。
    TDF_LabelSequence roots;
    shapeTool->GetFreeShapes(roots);
    std::vector<TopoDS_Shape> pieces;
    std::function<void(const TopoDS_Shape&)> collect = [&](const TopoDS_Shape& shape) {
        if (shape.IsNull()) return;
        if (shape.ShapeType() == TopAbs_COMPOUND) {
            // 默认迭代累计父位置与方向；统一容器自身保持单位变换。
            for (TopoDS_Iterator it(shape); it.More(); it.Next()) collect(it.Value());
        }
        else if (shape.ShapeType() == TopAbs_SHELL
            || shape.ShapeType() == TopAbs_SOLID
            || shape.ShapeType() == TopAbs_COMPSOLID) {
            // 仅收集完整壳/实体；不展开其内部面、边、顶点。
            pieces.push_back(shape);
        }
        // 独立 FACE/WIRE/EDGE/VERTEX 不参与修复和最终导出。
    };
    for (int i = 1; i <= roots.Length(); ++i) collect(shapeTool->GetShape(roots.Value(i)));
    std::stable_sort(pieces.begin(), pieces.end(), [](const TopoDS_Shape& a, const TopoDS_Shape& b) {
        return a.ShapeType() < b.ShapeType();
    });
    TopTools_MapOfShape covered;
    std::vector<TopoDS_Shape> uniquePieces;
    for (const auto& piece : pieces) {
        if (covered.Contains(piece)) continue;
        uniquePieces.push_back(piece);
        TopExp::MapShapes(piece, covered); // 包含在 solid/shell 中的重复根也不再加入。
    }
    if (uniquePieces.empty()) return -1;
    log << "Grouped geometry repair: freeRoots=" << roots.Length()
        << ", uniquePieces=" << uniquePieces.size() << std::endl;

    auto report = [&](const ShapeRepairPipeline::Result& result) {
        for (const auto& stage : result.stageReports) {
            log << "[Stage " << stage.name << "] accepted=" << stage.accepted
                << ", modified=" << stage.modified << ", repairCount=" << stage.repairCount
                << ", rejection reason=" << stage.rejectionReason << '\n';
            for (const auto& message : stage.messages) log << message << '\n';
        }
        for (const auto& message : result.messages) log << message << '\n';
        log << "Pipeline modified=" << result.modified << ", rolledBack=" << result.rolledBack
            << ", committed=" << result.totalRepairCount << std::endl;
    };
    BRep_Builder builder;
    TopoDS_Compound original, combined;
    builder.MakeCompound(original);
    builder.MakeCompound(combined);
    int total = 0;
    for (std::size_t i = 0; i < uniquePieces.size(); ++i) {
        const auto& piece = uniquePieces[i];
        builder.Add(original, piece);
        TopoDS_Shape accepted = piece;
        if (piece.ShapeType() <= TopAbs_SHELL) {
            log << "[Initial piece " << i << "] type="
                << TopologyUtils::ShapeTypeName(piece) << std::endl;
            std::cout << "Initial repair piece " << i + 1 << "/" << uniquePieces.size()
                << ": " << TopologyUtils::ShapeTypeName(piece) << std::endl;
            BRepBuilderAPI_Copy copy(piece, Standard_True, Standard_False);
            if (copy.Shape().IsNull()) return -1;
            const auto result = pipeline.perform(copy.Shape());
            report(result);
            if (result.modified && !result.rolledBack && !result.resultShape.IsNull()) {
                accepted = result.resultShape;
                total += result.totalRepairCount;
            }
        }
        builder.Add(combined, accepted);
    }
    ShapeRepairPipeline rounds;
    rounds.addStage(ShapeRepairStages::MakeConservativeSewShellsStage(0.26, true));

    TopoDS_Shape working = combined;
    for (int round = 1; round <= 1; ++round) {
        log << "[Global round " << round << "]" << std::endl;
        std::cout << "Global sewing/repair round " << round << "/1" << std::endl;
        BRepBuilderAPI_Copy copy(working, Standard_True, Standard_False);
        if (copy.Shape().IsNull()) return -1;
        const auto result = rounds.perform(copy.Shape());
        report(result);
        if (!result.modified || result.rolledBack || result.resultShape.IsNull()) {
            log << "Stop: no committed modification in round " << round << std::endl;
            break;
        }
        working = result.resultShape;
        total += result.totalRepairCount;
    }
    std::string error;
    const auto before = RepairDiagnostics::CollectBRepCheckFailures(original, &error);
    if (!error.empty()) { log << error; return -1; }
    const auto after = RepairDiagnostics::CollectBRepCheckFailures(working, &error);
    if (!error.empty()) { log << error; return -1; }
    const auto added = RepairDiagnostics::FindBRepCheckFailureCountIncrease(before, after);
    if (!added.empty() || TopologyUtils::CountFreeEdges(working) > TopologyUtils::CountFreeEdges(original)) {
        log << "Final validation failed; exporting original combined geometry. "
            << RepairDiagnostics::DescribeBRepCheckFailures(added) << std::endl;
        working = original;
        total = 0;
    }
    log << "Final committed operations=" << total << std::endl;
    for (const auto& entry : { std::make_pair("Original", TopoDS_Shape(original)),
                              std::make_pair("Final", working) }) {
        log << entry.first << " topology:";
        for (TopAbs_ShapeEnum type : { TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX }) {
            TopTools_IndexedMapOfShape shapes;
            TopExp::MapShapes(entry.second, type, shapes);
            log << ' ' << TopologyUtils::ShapeTypeName(type) << '=' << shapes.Extent();
        }
        log << " freeEdges=" << TopologyUtils::CountFreeEdges(entry.second) << std::endl;
    }
    std::vector<std::string> finalChecks;
    RepairDiagnostics::AppendShapeCheckSummary(finalChecks, "Final", working);
    for (const auto& message : finalChecks) log << message << std::endl;
    // 普通 STEP writer 不传入文档，且显式关闭装配导出。
    const int previousAssemblyMode = Interface_Static::IVal("write.step.assembly");
    Interface_Static::SetIVal("write.step.assembly", 0);
    const bool written = exportStepShape(working, filePath);
    Interface_Static::SetIVal("write.step.assembly", previousAssemblyMode);
    return written ? total : -1;
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
