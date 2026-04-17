#include "ApplicationInit.h"

#include <iostream>

#include <IGESCAFControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

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

    std::cout << "IGES file read successfully." << std::endl;
    return transferAndCacheTopLevelShapes(aReader.Transfer(aDoc_), "IGES");
}

bool ApplicationInit::exportStep(const TopoDS_Shape& shape, const std::string& filePath) const
{
    Interface_Static::SetIVal("write.step.assembly", 0);
    if (shape.IsNull()) {
        std::cerr << "Error: The shape to export is null." << std::endl;
        return false;
    }

    std::cout << "Writing STEP file to: " << filePath << "......" << std::endl;
    STEPControl_Writer writer;

    const IFSelect_ReturnStatus status = writer.Transfer(shape, STEPControl_AsIs);
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to transfer shape to STEP writer." << std::endl;
        return false;
    }

    if (writer.Write(filePath.c_str()) == IFSelect_RetDone) {
        std::cout << "STEP file written successfully." << std::endl;
        return true;
    }

    std::cerr << "Error: Failed to write STEP file." << std::endl;
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
    shapeTool->GetFreeShapes(freeShapes);

    TopTools_MapOfShape shapeMap;
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
