#include "ApplicationInit.h"

#include <iostream>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <TDF_ChildIterator.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <STEPControl_Writer.hxx>

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
    IFSelect_ReturnStatus status = aReader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to read STEP file." << std::endl;
        return false;
    }
    std::cout << "STEP file read successfully." << std::endl;
    if (aReader.Transfer(aDoc_)) {
        std::cout << "STEP file transferred to document successfully." << std::endl;
        return true;
    }
    else {
        std::cerr << "Error: Failed to transfer STEP file to document." << std::endl;
        return false;
    }
}


bool ApplicationInit::exportStep(const TopoDS_Shape& shape, const std::string& filePath) const
{
    if (shape.IsNull()) {
        std::cerr << "Error: The shape to export is null." << std::endl;
        return false;
    }

    std::cout << "Writing STEP file to: " << filePath << "......" << std::endl;
    STEPControl_Writer writer;

    // 将 Shape 转换入写入器
    IFSelect_ReturnStatus status = writer.Transfer(shape, STEPControl_AsIs);
    if (status != IFSelect_RetDone) {
        std::cerr << "Error: Failed to transfer shape to STEP writer." << std::endl;
        return false;
    }

    // 执行物理写入
    status = writer.Write(filePath.c_str());
    if (status == IFSelect_RetDone) {
        std::cout << "STEP file written successfully." << std::endl;
        return true;
    }
    else {
        std::cerr << "Error: Failed to write STEP file." << std::endl;
        return false;
    }
}

void ApplicationInit::traverseLabels(const TDF_Label& aLabel, std::vector<TopoDS_Shape>& shapeVec) const
{
    if (XCAFDoc_ShapeTool::IsShape(aLabel)) {
        TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aLabel);
        if (!aShape.IsNull()) {
            shapeVec.push_back(aShape);
        }
    }
    // 递归遍历子标签
    TDF_ChildIterator it(aLabel);
    for (; it.More(); it.Next()) {
        traverseLabels(it.Value(), shapeVec);
    }
}

std::vector<TopoDS_Shape> ApplicationInit::getAllShapes() const
{
    std::vector<TopoDS_Shape> shapeVec;
    if (aDoc_.IsNull()) return shapeVec;

    // 获取文档的 ShapeTool
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(aDoc_->Main());

    // 获取 Free Shapes (只有顶层装配体/独立的零件，不包含子层级)
    TDF_LabelSequence freeShapes;
    shapeTool->GetFreeShapes(freeShapes);

    // 遍历顶层标签并提取 Shape
    for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
        TopoDS_Shape aShape = shapeTool->GetShape(freeShapes.Value(i));
        if (!aShape.IsNull()) {
            shapeVec.push_back(aShape);
        }
    }
    return shapeVec;
}

std::vector<TopoDS_Solid> ApplicationInit::getAllSolids() const
{
    std::vector<TopoDS_Solid> solidVec;
    std::vector<TopoDS_Shape> allShapes = getAllShapes();

    for (const auto& shape : allShapes) {
        if (shape.ShapeType() == TopAbs_SOLID) {
            solidVec.push_back(TopoDS::Solid(shape));
        }
        else if (shape.ShapeType() == TopAbs_COMPOUND || shape.ShapeType() == TopAbs_COMPSOLID) {
            // 递归探索复合形状中的 Solid
            for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
                solidVec.push_back(TopoDS::Solid(exp.Current()));
            }
        }
    }
    return solidVec;
}

std::vector<TopoDS_Shell> ApplicationInit::getAllShells() const
{
    std::vector<TopoDS_Shell> shellVec;
    std::vector<TopoDS_Shape> allShapes = getAllShapes();

    for (const auto& shape : allShapes) {
        if (shape.ShapeType() == TopAbs_SHELL) {
            shellVec.push_back(TopoDS::Shell(shape));
        }
        else if (shape.ShapeType() == TopAbs_SOLID ||
            shape.ShapeType() == TopAbs_COMPOUND ||
            shape.ShapeType() == TopAbs_COMPSOLID) {
            // 递归探索形状中的 Shell
            for (TopExp_Explorer exp(shape, TopAbs_SHELL); exp.More(); exp.Next()) {
                shellVec.push_back(TopoDS::Shell(exp.Current()));
            }
        }
    }
    return shellVec;
}

std::vector<TopoDS_Edge> ApplicationInit::getAllEdges() const
{
    std::vector<TopoDS_Edge> edgeVec;
    std::vector<TopoDS_Shape> allShapes = getAllShapes();

    for (const auto& shape : allShapes) {
        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            edgeVec.push_back(TopoDS::Edge(exp.Current()));
        }
    }
    return edgeVec;
}