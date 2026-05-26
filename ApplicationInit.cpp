#include "ApplicationInit.h"

#include <cctype>
#include <iostream>
#include <IGESCAFControl_Reader.hxx>
#include <IGESData_BasicEditor.hxx>
#include <IGESData_GlobalSection.hxx>
#include <IGESData_IGESModel.hxx>
#include <Message.hxx>
#include <Message_Printer.hxx>
#include <Message_SequenceOfPrinters.hxx>
#include <STEPCAFControl_Writer.hxx>
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


namespace {
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

void ApplicationInit::collectTopLevelTargets(std::vector<EditableContainerTarget>& targets) const {
    if (aDoc_.IsNull()) {
        return;
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());
    if (shapeTool.IsNull()) {
        return;
    }
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

int ApplicationInit::healNarrowGapsOnLowestContainersAndExport(const std::string& filePath,
                                                               const HealNarrowGap::Params& params)
{
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
        const int containerRepairedCount = narrowGapHealer.Perform();

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

        for (const std::string& message : narrowGapHealer.getMessages()) {
            std::cout << "[HealNarrowGap] Container " << i << ": " << message << std::endl;
        }
    }

    shapeTool->UpdateAssemblies();
    refreshShapeCache();
    std::cout << "Replaced lowest containers after narrow-gap heal: " << replacedCount << std::endl;
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
