#include "ApplicationInit.h"

#include <iostream>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_ChildIterator.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>

ApplicationInit::ApplicationInit()
{
    anApp_ = new TDocStd_Application();
    anApp_->NewDocument("BinXCAF", aDoc_);
}

ApplicationInit::~ApplicationInit()
{
    // 文档由应用管理，这里不需要显式关闭
}

bool ApplicationInit::importStep(const std::string& filePath)
{
    STEPCAFControl_Reader aReader;
    std::cout << "Reading STEP file: " << filePath << "......" << std::endl;
    IFSelect_ReturnStatus status = aReader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone) {
        std::cout << "Error: Failed to read STEP file." << std::endl;
        return false;
    }
    std::cout << "STEP file read successfully." << std::endl;
    if (aReader.Transfer(aDoc_)) {
        std::cout << "STEP file transferred to document successfully." << std::endl;
        return true;
    }
    else {
        std::cout << "Error: Failed to transfer STEP file to document." << std::endl;
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
    if (!aDoc_.IsNull()) {
        traverseLabels(aDoc_->Main(), shapeVec);
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
