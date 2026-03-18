#pragma once
#include <string>
#include <vector>

#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Edge.hxx>
#include <TDF_Label.hxx>

class ApplicationInit {
public:
    ApplicationInit();
    ~ApplicationInit();

    // 获取 OCC 应用和文档
    Handle(TDocStd_Application) getApp() const { return anApp_; }
    Handle(TDocStd_Document) getDoc() const { return aDoc_; }

    // 导入 STEP 文件
    bool importStep(const std::string& filePath);

    // 获取所有 TopoDS_Shape（包括 Solid 和 Shell）
    std::vector<TopoDS_Shape> getAllShapes() const;

    // 获取所有 TopoDS_Solid
    std::vector<TopoDS_Solid> getAllSolids() const;

    // 获取所有 TopoDS_Shell
    std::vector<TopoDS_Shell> getAllShells() const;

    // 获取所有 TopoDS_Edge
    std::vector<TopoDS_Edge> getAllEdges() const;

private:
    void traverseLabels(const TDF_Label& aLabel, std::vector<TopoDS_Shape>& shapeVec) const;

    Handle(TDocStd_Application) anApp_;
    Handle(TDocStd_Document) aDoc_;
};
