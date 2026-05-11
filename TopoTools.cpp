#include "TopoTools.h"

#include <chrono>
#include <ctime>
#include <unordered_set>

#include <topo/body.h>
#include <topo/edge.h>
#include <topo/face.h>
#include <topo/vertex.h>

#include <TopExp.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomConvert.hxx>
#include <TopTools_HSequenceOfShape.hxx>

#include <Geom_Curve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include<Geom_ToroidalSurface.hxx>

#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <IGESControl_Reader.hxx>
#include <ShapeFix_Shape.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <IGESCAFControl_Reader.hxx>
#include <TDF_ChildIterator.hxx>

TopoTools::TopoTools(double linearDeflection, const std::string& fileName)
{
    aLinearDeflection_ = linearDeflection;
    std::string debugFilePath = std::format("./results/{}-debugLog.txt", fileName);
    fileName_ = fileName;
    aDebugLog_.open(debugFilePath);
    if (!aDebugLog_.is_open())
    {
        std::cout << "Failed to open debug log file." << std::endl;
    }
    auto now = std::chrono::system_clock::now();
    auto localTime = std::chrono::zoned_time{ std::chrono::current_zone(), now };
    aDebugLog_ << std::format("Start time: {:%Y-%m-%d %X}", localTime) << std::endl;

    // 初始化应用程序
    anApp_ = new TDocStd_Application();

    // 创建一个新文档
    anApp_->NewDocument("BinXCAF", aDoc_);


}

TopoTools::~TopoTools()
{
    aDebugLog_ << "Done!" << std::endl;
    aDebugLog_.close();
    // new 出来的数据后面记得删除
    // TODO: 看实际情况是否使用智能指针，一个是OCC自带的Handle，一个是C++的智能指针
}


// 将一个TopoDS_Solid转换为Body
// 此方法会导致收录重复的顶点，边因为是指针所以没重复的
Body* TopoTools::convertToBody(const TopoDS_Solid& aTopoSolid)
{
    Body* aBody = new Body();

    NCollection_DataMap<TopoDS_Edge, Edge*> aTopoEdge_EdgeMap;
    NCollection_DataMap<TopoDS_Wire, Loop*> aTopoWire_LoopMap;

    aDebugLog_ << "Converting a new body......" << std::endl;
    // 遍历TopoDS_Solid中的所有边，将TopoDS_Edge转换为Edge，并保存到aTopoEdge_EdgeMap中
    for (TopExp_Explorer anEdgeExp(aTopoSolid, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
    {
        TopoDS_Edge aTopoEdge = TopoDS::Edge(anEdgeExp.Current());
        Edge* anEdge = convertToEdge(aTopoEdge);
        aTopoEdge_EdgeMap.Bind(aTopoEdge, anEdge);
    }
    // 遍历TopoDS_Solid中的所有环，并保存到aTopoWire_LoopMap中
    for (TopExp_Explorer aWireExp(aTopoSolid, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
    {
        TopoDS_Wire aTopoWire = TopoDS::Wire(aWireExp.Current());
        Loop* aLoop = new Loop();
        for (TopExp_Explorer aExp(aTopoWire, TopAbs_EDGE); aExp.More(); aExp.Next())
        {
            TopoDS_Edge aTopoEdge = TopoDS::Edge(aExp.Current());
            Edge* anEdge = aTopoEdge_EdgeMap.Find(aTopoEdge);
            aLoop->addEdge(*anEdge);
        }
        aTopoWire_LoopMap.Bind(aTopoWire, aLoop);
    }
    // 遍历TopoDS_Solid中的所有面，并保存到aBody中
    for (TopExp_Explorer aFaceExp(aTopoSolid, TopAbs_FACE); aFaceExp.More(); aFaceExp.Next(), faceIndex_++)
    {
        TopoDS_Face aTopoFace = TopoDS::Face(aFaceExp.Current());
        Face* aFace = convertToFace(aTopoFace);
        aFace->index = faceIndex_;
        for (TopExp_Explorer aWireExp(aTopoFace, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
        {
            TopoDS_Wire aTopoWire = TopoDS::Wire(aWireExp.Current());
            Loop* aLoop = aTopoWire_LoopMap.Find(aTopoWire);
            aFace->addLoop(*aLoop);
        }
        for (TopExp_Explorer anEdgeExp(aTopoFace, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
        {
            TopoDS_Edge aTopoEdge = TopoDS::Edge(anEdgeExp.Current());
            Edge* anEdge = aTopoEdge_EdgeMap.Find(aTopoEdge);
            aFace->addEdge(*anEdge);
        }
        aBody->addFace(*aFace);
    }

    return aBody;
}

Body* TopoTools::convertToBody(const std::vector<TopoDS_Shell>& aTopoShellVec)
{
    Body* aBody = new Body();

    // 模版推导有问题
    //std::map<TopoDS_Edge, Edge*> aTopoEdge_EdgeMap;
    //std::map<TopoDS_Wire, Loop*> aTopoWire_LoopMap;

    NCollection_DataMap<TopoDS_Edge, Edge*> aTopoEdge_EdgeMap;
    NCollection_DataMap<TopoDS_Wire, Loop*> aTopoWire_LoopMap;


    std::unordered_set<TopoDS_Edge> aTopoEdgeSet;
    std::unordered_set<TopoDS_Wire> aTopoWireSet;
    std::unordered_set<TopoDS_Face> aTopoFaceSet;

    std::vector<TopoDS_Edge> aTopoEdgeVec;
    std::vector<TopoDS_Wire> aTopoWireVec;
    std::vector<TopoDS_Face> aTopoFaceVec;

    aDebugLog_ << "Converting a new body......" << std::endl;
    // 遍历aTopoShellVec中的所有边，并保存到aTopoEdge_EdgeMap中
    for (auto& aTopoShell : aTopoShellVec)
    {
        for (TopExp_Explorer anEdgeExp(aTopoShell, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
        {
            TopoDS_Edge aTopoEdge = TopoDS::Edge(anEdgeExp.Current());
            if (aTopoEdgeSet.find(aTopoEdge) == aTopoEdgeSet.end()) // 如果集合中没这个Edge
            {
                aTopoEdgeSet.insert(aTopoEdge);
                Edge* anEdge = convertToEdge(aTopoEdge);
                aTopoEdge_EdgeMap.Bind(aTopoEdge, anEdge);
            }
        }
    }

    // 遍历aTopoShellVec中的所有环，并保存到aTopoWire_LoopMap中
    for (auto& aTopoShell : aTopoShellVec)
    {
        for (TopExp_Explorer aWireExp(aTopoShell, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
        {
            TopoDS_Wire aTopoWire = TopoDS::Wire(aWireExp.Current());
            if (aTopoWireSet.find(aTopoWire) == aTopoWireSet.end())
            {
                aTopoWireSet.insert(aTopoWire);
                Loop* aLoop = new Loop();
                for (TopExp_Explorer aExp(aTopoWire, TopAbs_EDGE); aExp.More(); aExp.Next())
                {
                    TopoDS_Edge aTopoEdge = TopoDS::Edge(aExp.Current());
                    Edge* anEdge = aTopoEdge_EdgeMap.Find(aTopoEdge);
                    aLoop->addEdge(*anEdge);
                }
                aTopoWire_LoopMap.Bind(aTopoWire, aLoop);
            }
        }
    }
    // 遍历aTopoShellVec中的所有面，并保存到aBody中
    for (auto& aTopoShell : aTopoShellVec)
    {
        for (TopExp_Explorer aFaceExp(aTopoShell, TopAbs_FACE); aFaceExp.More(); aFaceExp.Next(), faceIndex_++)
        {
            TopoDS_Face aTopoFace = TopoDS::Face(aFaceExp.Current());
            if (aTopoFaceSet.find(aTopoFace) == aTopoFaceSet.end())
            {
                aTopoFaceSet.insert(aTopoFace);
                // Face不会重复
                // 那么什么会重复？
                Face* aFace = convertToFace(aTopoFace);
                aFace->index = faceIndex_;
                for (TopExp_Explorer aWireExp(aTopoFace, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
                {
                    TopoDS_Wire aTopoWire = TopoDS::Wire(aWireExp.Current());
                    Loop* aLoop = aTopoWire_LoopMap.Find(aTopoWire);
                    aFace->addLoop(*aLoop);
                }
                for (TopExp_Explorer anEdgeExp(aTopoFace, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
                {
                    TopoDS_Edge aTopoEdge = TopoDS::Edge(anEdgeExp.Current());
                    Edge* anEdge = aTopoEdge_EdgeMap.Find(aTopoEdge);
                    aFace->addEdge(*anEdge);
                }
                aBody->addFace(*aFace);
            }
        }
    }

    aBody->isSolid = false;
    return aBody;
}



Vec3d TopoTools::convertToVec3d(const TopoDS_Vertex& aOCC_TopoVertex)
{
    gp_Pnt aGpPoint = BRep_Tool::Pnt(aOCC_TopoVertex);
    Vec3d aVec3d(aGpPoint.X(), aGpPoint.Y(), aGpPoint.Z());
    return aVec3d;
}


Vertex TopoTools::convertToVertex(const TopoDS_Vertex& aOCC_TopoVertex)
{
    gp_Pnt aGpPoint = BRep_Tool::Pnt(aOCC_TopoVertex);
    Vertex aVertex(aGpPoint.X(), aGpPoint.Y(), aGpPoint.Z());
    return aVertex;
}


// 此函数只填充了Edge类中的curve属性，其他属性需要手动填充
Edge* TopoTools::convertToEdge(const TopoDS_Edge& anOCC_TopoEdge)
{
    Edge* anEdge = new Edge();

    TopoDS_Vertex aTopoVertex1 = TopExp::FirstVertex(anOCC_TopoEdge);
    TopoDS_Vertex aTopoVertex2 = TopExp::LastVertex(anOCC_TopoEdge);

    anEdge->vertex1 = convertToVertex(aTopoVertex1);
    anEdge->vertex2 = convertToVertex(aTopoVertex2);

    anEdge->curve = convertToNurbsCurve(anOCC_TopoEdge);

    return anEdge;
}


// 此函数只填充了Face类中的geom_surface属性和discrete_surface属性，其他属性需要手动填充
Face* TopoTools::convertToFace(const TopoDS_Face& anOCC_TopoFace)
{
    BRepMesh_IncrementalMesh aMesh(anOCC_TopoFace, aLinearDeflection_);
    Face* aFace = new Face();
    aFace->geom_surface = convertToNurbsSurface(anOCC_TopoFace);
    aFace->discrete_surface = convertToDiscreteSurface(anOCC_TopoFace);
    return aFace;
}


// 提取TopoDS_Edge的几何曲线，将其转换为NurbsCurve类型，并将其保存到自定义的Edge类中
NurbsCurve* TopoTools::convertToNurbsCurve(const TopoDS_Edge& anOCC_TopoEdge)
{
    // 将TopoDS_Edge类型转换为Geom_Curve类型，并转换为Geom_BSplineCurve类型
    double aFirstParam = 0, aLastParam = 0;
    TopLoc_Location aLocation;
    Handle(Geom_Curve) aGeomCurve = BRep_Tool::Curve(anOCC_TopoEdge, aLocation, aFirstParam, aLastParam);
    if (aGeomCurve.IsNull())
    {
        if (BRep_Tool::Degenerated(anOCC_TopoEdge))
        {
            aDebugLog_ << "Warning: This edge is degenerated and has no 3D curve." << std::endl;
        }
        else
        {
            aDebugLog_ << "Error: Failed to retrieve 3D curve for the edge (possibly missing geometry)." << std::endl;
        }
        return nullptr;
    }
    // 打印提取出的曲线的 类型 和 参数范围
    aDebugLog_ << "Curve type: " << aGeomCurve->DynamicType()->Name() << std::endl;
    //aDebugLog_ << "FirstParam: " << aFirstParam << ", LastParam: " << aLastParam << std::endl;
    Handle(Geom_BSplineCurve) anOCC_NurbsCurve;
    // 如果为直线，则需要对他范围进行限定
    if (!aGeomCurve->IsKind("Geom_BSplineCurve"))
    {
        // 根据BRep_Tool获得参数范围，将该直线裁剪，那么会得到Geom_TrimmedCurve类型的对象
        aGeomCurve = new Geom_TrimmedCurve(aGeomCurve, aFirstParam, aLastParam);
    }

    anOCC_NurbsCurve = GeomConvert::CurveToBSplineCurve(aGeomCurve);


    // 创建nubrs曲线的属性
    std::vector<Vec3d> controlPoints;
    std::vector<double> weights{ 0 };
    int order;
    std::vector<double> knots{ 0 };

    // 获取Geom_BSplineCurve的属性，并保存到自定义的Edge类中
    order = anOCC_NurbsCurve->Degree() + 1;
    // 存放控制点和权重
    //aDebugLog_ << "Poles and Weights: " <<std::endl;
    for (int i = 1; i <= anOCC_NurbsCurve->NbPoles(); i++)
    {
        gp_Pnt aGpPoint = anOCC_NurbsCurve->Pole(i);
        controlPoints.push_back(Vec3d(aGpPoint.X(), aGpPoint.Y(), aGpPoint.Z()));
        weights.push_back(anOCC_NurbsCurve->Weight(i));
        //aDebugLog_ << aGpPoint.X() << "\t " << aGpPoint.Y() << "\t " << aGpPoint.Z() << "\t Weight: " << nurbsCurve->Weight(i) << std::endl;
    }

    // 存放节点
    TColStd_Array1OfReal aKnotArray = anOCC_NurbsCurve->KnotSequence();
    //aDebugLog_ << "KnotSequence: " << std::endl;
    for (int i = aKnotArray.Lower(); i <= aKnotArray.Upper(); i++)
    {
        //aDebugLog_ << aKnotArray(i) << "\t ";
        knots.push_back(aKnotArray(i));
    }
    NurbsCurve* aNurbsCurve = new NurbsCurve(controlPoints, weights, order, knots);
    return aNurbsCurve;
}


// 提取TopoDS_Face几何曲线，将其转换为NurbsSurface类型，并将其保存到自定义的Face类中
NurbsSurface* TopoTools::convertToNurbsSurface(const TopoDS_Face& anOCC_TopoFace)
{
    // 将TopoDS_Face类型转换为Geom_Surface类型，并转换为Geom_BSplineSurface类型
    TopLoc_Location aLocation;
    Handle(Geom_Surface) anOCC_GeomSurface = BRep_Tool::Surface(anOCC_TopoFace, aLocation);
    if (anOCC_GeomSurface.IsNull())
    {
        aDebugLog_ << "Error: No geometric surface found for the face." << std::endl;
        return nullptr;
    }
    if (!aLocation.IsIdentity()) 
    {
        anOCC_GeomSurface = Handle(Geom_Surface)::DownCast(anOCC_GeomSurface->Transformed(aLocation.Transformation()));
    }

	std::string surfaceType = anOCC_GeomSurface->DynamicType()->Name();
    aDebugLog_ << "Surface type: " << surfaceType << std::endl;
    Handle(Geom_BSplineSurface) anOCC_BsplineSurface;
    if (anOCC_GeomSurface->IsKind(STANDARD_TYPE(Geom_BSplineSurface))) {
        anOCC_BsplineSurface = Handle(Geom_BSplineSurface)::DownCast(anOCC_GeomSurface);
    }
    else if (anOCC_GeomSurface->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface))) {
        Handle(Geom_RectangularTrimmedSurface) anOCC_TrimmmedSurface = Handle(Geom_RectangularTrimmedSurface)::DownCast(anOCC_GeomSurface);
		Handle(Geom_Surface) aBaseSurface = anOCC_TrimmmedSurface->BasisSurface();
        if (aBaseSurface->IsKind(STANDARD_TYPE(Geom_ToroidalSurface))) {
            anOCC_BsplineSurface = GeomConvert::SurfaceToBSplineSurface(aBaseSurface);
        }
        else {
            anOCC_BsplineSurface = GeomConvert::SurfaceToBSplineSurface(anOCC_TrimmmedSurface);
        }
    }
    else{
        //aDebugLog_ << "Warning: It's not a B-spline surface." << std::endl;
        double umin, umax, vmin, vmax;
        double surfUmin, surfUmax, surfVmin, surfVmax;
        const double aMinRange = 1e-6;  // 容差值，根据需要调整

        anOCC_GeomSurface->Bounds(surfUmin, surfUmax, surfVmin, surfVmax);
        BRepTools::UVBounds(anOCC_TopoFace, umin, umax, vmin, vmax);

        // TODO: 考虑处理边界情况，但这种可能欠妥
        umin = Max(umin, surfUmin);
        umax = Min(umax, surfUmax);
        vmin = Max(vmin, surfVmin);
        vmax = Min(vmax, surfVmax);

        // TODO: 有时候会出现超过2π的情况，导致程序报错
        // 可能原因: 浮点数精度问题造成的超过2π
        if (Abs(umax - umin) < aMinRange || Abs(vmax - vmin) < aMinRange)
        {
            aDebugLog_ << "Error: the surface is too small" << std::endl;
            return nullptr;
        }
        if (umin > umax || vmin > vmax) {
            aDebugLog_ << "Error: u,v parameters' range is wrong" << std::endl;
            return nullptr;
        }
        Handle(Geom_RectangularTrimmedSurface) anOCC_TrimmedSurface;
        anOCC_TrimmedSurface = new Geom_RectangularTrimmedSurface(anOCC_GeomSurface, umin, umax, vmin, vmax);
        anOCC_BsplineSurface = GeomConvert::SurfaceToBSplineSurface(anOCC_TrimmedSurface);
        if (anOCC_BsplineSurface.IsNull()) 
        {
            aDebugLog_ << "Error: can't convert this surface to Nurbs surface" << std::endl;
            return nullptr;
        }
    }

    std::vector<std::vector<Vec3d>> controlPoints;
    std::vector<double> weights{ 0 };
    int uOrder, vOrder;
    std::vector<double> uKnots{ 0 }, vKnots{ 0 };

    // 获取Geom_BSplineSurface的属性，并保存到自定义的Face类中
    uOrder = anOCC_BsplineSurface->UDegree() + 1;
    vOrder = anOCC_BsplineSurface->VDegree() + 1;
    for (int i = 1; i <= anOCC_BsplineSurface->NbUPoles(); i++)
    {
        std::vector<Vec3d> aRow;
        for (int j = 1; j <= anOCC_BsplineSurface->NbVPoles(); j++)
        {
            gp_Pnt aGpPoint = anOCC_BsplineSurface->Pole(i, j);
            aRow.push_back(Vec3d(aGpPoint.X(), aGpPoint.Y(), aGpPoint.Z()));
            weights.push_back(anOCC_BsplineSurface->Weight(i, j));
        }
        controlPoints.push_back(aRow);
    }

    // 获取节点
    TColStd_Array1OfReal aUKnotArray = anOCC_BsplineSurface->UKnotSequence();
    for (int i = aUKnotArray.Lower(); i <= aUKnotArray.Upper(); i++)
    {
        uKnots.push_back(aUKnotArray(i));
    }

    TColStd_Array1OfReal aVKnotArray = anOCC_BsplineSurface->VKnotSequence();
    for (int i = aVKnotArray.Lower(); i <= aVKnotArray.Upper(); i++)
    {
        vKnots.push_back(aVKnotArray(i));
    }

    NurbsSurface* aNurbsSurface = new NurbsSurface(controlPoints, weights, uOrder, vOrder, uKnots, vKnots);
    return aNurbsSurface;
}





// 提取离散TopoDS_Face的离散数据，将其转换为DiscreteSurface类型，并将其保存到自定义的DiscreteSurface类中
DiscreteSurface* TopoTools::convertToDiscreteSurface(const TopoDS_Face& anOCC_TopoFace)
{
    int elementsNum = 0;
    int pointsNum = 0;
    std::vector<Point> aPointVec;
    std::vector<int> aPointIndexVec;
    std::vector<int> aElementPointsNumVec;

    TopLoc_Location aLocation;
    Handle(Poly_Triangulation) aTriangulation = BRep_Tool::Triangulation(anOCC_TopoFace, aLocation);
    if (!aTriangulation.IsNull())
    {
        // 获取顶点和三角形索引
        int numNodes = aTriangulation->NbNodes();
        pointsNum = numNodes;
        elementsNum = aTriangulation->NbTriangles();

        // 处理顶点和三角形
        // Node∈[1, numNodes]
        for (int i = 1; i <= numNodes; ++i)
        {
            gp_Pnt aGpPoint = aTriangulation->Node(i);
            if (!aLocation.IsIdentity()) {
                // 只有非恒等变换时才执行变换
                aGpPoint.Transform(aLocation.Transformation());
            }
            Point aPoint(aGpPoint.X(), aGpPoint.Y(), aGpPoint.Z());
            aPointVec.push_back(aPoint);
        }

        for (int i = 1; i <= aTriangulation->NbTriangles(); ++i)
        {
            Poly_Triangle aTriangle = aTriangulation->Triangle(i);
            int p1_index, p2_index, p3_index;
            aTriangle.Get(p1_index, p2_index, p3_index);
            // 处理三角形 (p1_index, p2_index, p3_index)
            aPointIndexVec.push_back(p1_index-1);
            aPointIndexVec.push_back(p2_index-1);
            aPointIndexVec.push_back(p3_index-1);
            aElementPointsNumVec.push_back(3);
        }
    }
    else
    {
        aDebugLog_ << "Error: No triangulation found for the face." << std::endl;
    }
    DiscreteSurface* aDiscreteSurface =
        new DiscreteSurface(elementsNum, pointsNum, aPointVec, aPointIndexVec, aElementPointsNumVec);

    return aDiscreteSurface;
}

/// @brief 根据输入的标签，从文件中获取对应的实体
/// @param aLabel 传入参数
/// @param aSolidVec 传出参数
/// @param aShellVec 传出参数
/// @param aShellLabelVec 传出参数
void TopoTools::traverseLabels(const TDF_Label& aLabel, std::vector<TopoDS_Solid>& aSolidVec, std::vector<TopoDS_Shell>& aShellVec, std::vector<TDF_Label>& aShellLabelVec)
{
    if (XCAFDoc_ShapeTool::IsShape(aLabel))
    {
        TopoDS_Shape aTopoShape = XCAFDoc_ShapeTool::GetShape(aLabel);
        if (aTopoShape.ShapeType() == TopAbs_SOLID)
        {

            aSolidVec.push_back(TopoDS::Solid(aTopoShape));
        }
        else if (aTopoShape.ShapeType() == TopAbs_SHELL)
        {
            aShellVec.push_back(TopoDS::Shell(aTopoShape));
            aShellLabelVec.push_back(aLabel);
        }

    }
    // 递归遍历当前标签的所有子标签
    TDF_ChildIterator it(aLabel);
    for (; it.More(); it.Next()) {
        traverseLabels(it.Value(), aSolidVec, aShellVec, aShellLabelVec);
    }
}


// 从STEP文件中导入shapes
bool TopoTools::importSTEP(std::string_view aFilePath)
{
    std::cout << "Reading STEP file: " << aFilePath << "......" << std::endl;
    STEPCAFControl_Reader aReader_STEP;
    auto t0 = std::chrono::steady_clock::now();
    const IFSelect_ReturnStatus aStatus = aReader_STEP.ReadFile(aFilePath.data());
    auto t1 = std::chrono::steady_clock::now();
    if (aStatus != IFSelect_RetDone)
    {
        std::cout << "A BREP file was not read successfully" << std::endl;
        return false;
    }
    std::cout << "Reading STEP file successfully." << std::endl;
    std::cout << "ReadFile cost: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << " ms" << std::endl;
    aReader_STEP.Transfer(aDoc_);

    auto t2 = std::chrono::steady_clock::now();
    std::cout << "Transfer cost: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
        << " ms" << std::endl;
    std::cout << "Transfer to doc successfully." << std::endl;
    return true;
}


bool TopoTools::importIGES(std::string_view aFilePath)
{
    std::cout << "Reading IGES file......" << aFilePath << std::endl;
    IGESCAFControl_Reader aReader_IGES;
    const IFSelect_ReturnStatus aStatus = aReader_IGES.ReadFile(aFilePath.data());
    if (aStatus != IFSelect_RetDone)
    {
        std::cout << "A BREP file was not read successfully" << std::endl;
        return false;
    }
    aReader_IGES.Transfer(aDoc_);
    return true;
}


/// @brief 对外的接口，获取所有转换后的 Body 对象，包括由实体和壳体生成的体。
/// @return 一个包含所有 Body 指针的 std::vector。
std::vector<Body*> TopoTools::getBodies()
{
    std::cout << "Converting solids to bodies......" << std::endl;
    std::vector<TopoDS_Solid> aSolidVec;
    std::vector<TopoDS_Shell> aShellVec;
    std::vector<TDF_Label> aShellLabelVec;
    traverseLabels(aDoc_->Main(), aSolidVec, aShellVec, aShellLabelVec);
	std::cout << "Finished traversing labels." << std::endl;
    int aSolidIndex = 1;
    for (auto& aTopoShape : aSolidVec)
    {
        aDebugLog_ << std::format("\nSolid[{}]:", aSolidIndex) << std::endl;
        TopoDS_Solid aTopoSolid = TopoDS::Solid(aTopoShape);
        Body* aBody = convertToBody(aTopoSolid);
        aBodyVec_.push_back(aBody);
        aSolidIndex++;
    }

    // 思路：
    // 将这个aShellSeq进行遍历，然后进行分组，有相同父亲的分为一组
    std::vector<std::vector<TopoDS_Shell>> aSortedShellVec;
    std::vector<bool> isSelected(aShellVec.size(), false);
    for (int i = 0; i < aShellVec.size(); i++)
    {
        if (isSelected[i] == true)
            continue;
        else
        {
            std::vector<TopoDS_Shell> aTempVec;
            aTempVec.push_back(TopoDS::Shell(aShellVec[i]));
            isSelected[i] = true;
            for (int j = i + 1; j < aShellVec.size(); j++)
            {
                // 如果上级标签一致，代表父亲一致
                if (aShellLabelVec[j].Father() == aShellLabelVec[i].Father())
                {
                    isSelected[j] = true;
                    aTempVec.push_back(TopoDS::Shell(aShellVec[j]));
                }
            }
            aSortedShellVec.push_back(aTempVec);
        }
    }


    aDebugLog_ << "The following solids is generated by shells" << std::endl << std::endl;
    for (auto& aShellVec : aSortedShellVec)
    {
        aDebugLog_ << std::format("\nSolid[{}]:", aSolidIndex) << std::endl;
        Body* aBody = convertToBody(aShellVec);
        aBodyVec_.push_back(aBody);
        aSolidIndex++;
    }

    if (aBodyVec_.size() == 0)
    {
        std::cout << "No body found" << std::endl;
        return aBodyVec_;
    }
    //else
    //{
    //    std::cout << std::format("Solid number in the compound: {}", aBodyVec_.size()) << std::endl;
    //}

    return aBodyVec_;
}


bool TopoTools::exportVTK(std::string_view aFileName, const std::vector<Face*>& aFaceVec)
{
    std::cout << "Exporting to VTK file......" << std::endl;
    std::string aFileNameWithPath = std::format("{}.vtk", aFileName.data());
    std::ofstream aFileStream(aFileNameWithPath);

    if (!aFileStream.is_open())
        return false;

    aFileStream << "# vtk DataFile Version 3.0" << std::endl; // 此为文件版本和标识符
    aFileStream << "background mesh" << std::endl; // 此为标题
    aFileStream << "ASCII" << std::endl; // 文件格式
    aFileStream << "DATASET UNSTRUCTURED_GRID" << std::endl; // 数据集格式，此为 非结构化网格

    int elementNum = 0;
    int pointNum = 0;

    for (auto& aFace : aFaceVec)
    {
        auto& aDiscreteSurface = aFace->discrete_surface;
        elementNum += aDiscreteSurface->get_elements_num();
        pointNum += aDiscreteSurface->get_points_num();
    }
    aFileStream << std::format("POINTS {} double", pointNum) << std::endl; // 输出顶点数量

    // 输出顶点坐标
    for (auto& aFace : aFaceVec)
    {
        auto& aDiscreteSurface = aFace->discrete_surface;
        int tempPointNum = aDiscreteSurface->get_points_num();
        for (int i = 0; i < tempPointNum; i++)
        {
            Point aTempPoint = aDiscreteSurface->get_point(i);
            aFileStream << std::format("{:.5f} {:.5f} {:.5f}", aTempPoint[0], aTempPoint[1], aTempPoint[2]) << std::endl;
        }
    }

    // 输出一个单元的数量
    aFileStream << std::format("CELLS {} {}", elementNum, elementNum * 4) << std::endl;
    int count = 0;
    for (auto& aFace : aFaceVec)
    {
        auto& aDiscreteSurface = aFace->discrete_surface;
        int tempElementNum = aDiscreteSurface->get_elements_num();
        for (int i = 0; i < tempElementNum * 3; i += 3)
        {
            // 在这里一个单元为三角形，所以需要输出三个顶点索引
            aFileStream << std::format("3 {} {} {}",
                count + aDiscreteSurface->get_point_index(i),
                count + aDiscreteSurface->get_point_index(i + 1),
                count + aDiscreteSurface->get_point_index(i + 2)
            );
            aFileStream << std::endl;
        }
        count += aDiscreteSurface->get_points_num();
    }

    aFileStream << std::format("CELL_TYPES {} ", elementNum) << std::endl;
    for (int i = 0; i < elementNum; i++) {
        aFileStream << std::format("{}", 5) << std::endl;
    }

    // 数据集属性格式
    aFileStream << std::format("CELL_DATA {}", elementNum) << std::endl;
    aFileStream << std::format("SCALARS surfaceid int 1") << std::endl;
    aFileStream << std::format("LOOKUP_TABLE default") << std::endl;

    for (int i = 0; i < aFaceVec.size(); i++)
    {
        int tempElementNum = aFaceVec[i]->discrete_surface->get_elements_num();
        for (int j = 0; j < tempElementNum; j++)
            aFileStream << i << std::endl;
    }
    aFileStream.close();

    return true;
}


bool TopoTools::exportOBJ(std::string_view aFileName, const std::vector<Face*>& aFaceVec)
{
    std::cout << "Exporting to OBJ file......" << std::endl;
    std::string aFileNameWithPath = std::format("{}.obj", aFileName.data());
    std::ofstream aFileStream(aFileNameWithPath);

    if (!aFileStream.is_open())
        return false;

    // 写入顶点
    for (auto& aFace : aFaceVec)
    {
        auto& aDiscreteSurface = aFace->discrete_surface;
        for (int i = 0; i < aDiscreteSurface->get_points_num(); i++)
        {
            Point aPoint = aDiscreteSurface->get_point(i);
            aFileStream << std::format("v {} {} {}",
                aPoint[0], aPoint[1], aPoint[2]) << std::endl;
        }
    }

    // 写入面
    int count = 1; // OBJ索引从1开始，这里只输出了顶点索引
    for (auto& aFace : aFaceVec)
    {
        auto& aDiscreteSurface = aFace->discrete_surface;
        int tempElementNum = aDiscreteSurface->get_elements_num();
        for (int i = 0; i < tempElementNum * 3; i += 3)
        {
            // 在这里一个单元为三角形，所以需要输出三个顶点索引
            aFileStream << std::format("f {} {} {}",
                count + aDiscreteSurface->get_point_index(i),
                count + aDiscreteSurface->get_point_index(i + 1),
                count + aDiscreteSurface->get_point_index(i + 2)
            );
            aFileStream << std::endl;
        }
        count += aDiscreteSurface->get_points_num();
    }

    return true;
}

