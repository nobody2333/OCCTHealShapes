#pragma once

#include <vector>
#include <format>
#include <string_view>

#include <geom/DiscreteCurve.h>
#include <geom/DiscreteSurface.h>
#include <geom/NurbsCurve.h>
#include <geom/NurbsSurface.h>

#include <topo/body.h>
#include <topo/edge.h>
#include <topo/face.h>
#include <topo/vertex.h>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Compound.hxx>

#include <NCollection_DataMap.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDocStd_Application.hxx>

class TopoTools // 该类名字有待商榷
{
public:
    TopoTools(double linearDeflection, const string& fileName);
    ~TopoTools();

    std::vector<Body*> getBodies();

    bool importSTEP(std::string_view aFilePath);
    bool importIGES(std::string_view aFilePath);

    bool exportVTK(std::string_view aFileName, const std::vector<Face*>& aFaceVec);
    bool exportOBJ(std::string_view aFileName, const std::vector<Face*>& aFaceVec);

private:
    Body* convertToBody(const TopoDS_Solid& anOCC_TopoSolid);
    Body* convertToBody(const std::vector<TopoDS_Shell>& anOCC_TopoShellVec);
    Face* convertToFace(const TopoDS_Face& anOCC_TopoFace);
    Edge* convertToEdge(const TopoDS_Edge& anOCC_TopoEdge);
    Vertex convertToVertex(const TopoDS_Vertex& anOCC_TopoVertex);
    Vec3d convertToVec3d(const TopoDS_Vertex& anOCC_TopoVertex);

    NurbsCurve* convertToNurbsCurve(const TopoDS_Edge& anOCC_TopoEdge);
    DiscreteCurve* convertToDiscreteCurve(const TopoDS_Edge& anOCC_TopoEdge);

    NurbsSurface* convertToNurbsSurface(const TopoDS_Face& anOCC_TopoFace);
    DiscreteSurface* convertToDiscreteSurface(const TopoDS_Face& anOCC_TopoFace);
    
    void traverseLabels(const TDF_Label& aLabel, std::vector<TopoDS_Solid>& aSolidVec, std::vector<TopoDS_Shell>& aShellVec, std::vector<TDF_Label>& aShellLabelVec);

private:
    double aLinearDeflection_;         // 线性偏差，控制着形状表面的最大允许偏离。在生成网格的过程中，任何从理论几何形状到实际网格表示之间的距离都不会超过这个值。
    int faceIndex_ = 1;
    string fileName_;
    Handle(TDocStd_Application) anApp_;
    Handle(XCAFDoc_ShapeTool) aShapeTool_;
    Handle(TDocStd_Document) aDoc_;
    
    std::vector<Body*> aBodyVec_;
    std::ofstream aDebugLog_;

};