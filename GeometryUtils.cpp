#include "GeometryUtils.h"

#include <BRep_Tool.hxx>
#include <BRepBndLib.hxx>
#include <TopExp.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <Bnd_Box.hxx>

#include <stdexcept>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double GeometryUtils::ComputeEdgeLength(const TopoDS_Edge& edge)
{
    if (edge.IsNull()) {
        return 0.0;
    }

    try {
        BRepAdaptor_Curve curveAdaptor(edge);
        return GCPnts_AbscissaPoint::Length(curveAdaptor, curveAdaptor.FirstParameter(), curveAdaptor.LastParameter());
    }
    catch (Standard_Failure&) {
        // 退化边或无几何边，回退到两点距离
        TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
        TopoDS_Vertex v2 = TopExp::LastVertex(edge);
        if (v1.IsNull() || v2.IsNull()) 
            return 0.0;
        return BRep_Tool::Pnt(v1).Distance(BRep_Tool::Pnt(v2));
    }
}

double GeometryUtils::ComputeIntersectionAngle(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex commonV)
{
    if (e1.IsNull() || e2.IsNull() || commonV.IsNull()) 
        return -1.0;

    try {
        // 使用 BRepAdaptor_Curve 会自动处理 TopoDS_Edge的Orientation和Location
        BRepAdaptor_Curve bac1(e1);
        BRepAdaptor_Curve bac2(e2);

        gp_Pnt pnt1, pnt2;
        gp_Vec d1, d2;

        // 判断公共顶点在 e1 上的位置，并获取指向外部的切向量
        if (commonV.IsSame(TopExp::FirstVertex(e1))) {
            bac1.D1(bac1.FirstParameter(), pnt1, d1);
            d1.Reverse(); // 在起点，导数方向指向内部，反转使其向外
        }
        else if (commonV.IsSame(TopExp::LastVertex(e1))) {
            bac1.D1(bac1.LastParameter(), pnt1, d1);
            // 在终点，导数方向顺着曲线向外，无需反转
        }
        else {
            return -1.0;
        }

        // 判断公共顶点在 e2 上的位置，并获取指向外部的切向量
        if (commonV.IsSame(TopExp::FirstVertex(e2))) {
            bac2.D1(bac2.FirstParameter(), pnt2, d2);
            d2.Reverse();
        }
        else if (commonV.IsSame(TopExp::LastVertex(e2))) {
            bac2.D1(bac2.LastParameter(), pnt2, d2);
        }
        else {
            return -1.0;
        }

        // 检查向量有效性
        if (d1.Magnitude() <= Precision::Confusion() || d2.Magnitude() <= Precision::Confusion()) {
            return -1.0;
        }

        return d1.Angle(d2);
    }
    catch (Standard_Failure&) {
        return -1.0;
    }
}

double GeometryUtils::ComputeBoundingBoxDiagonal(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) {
        return 0.0;
    }

    Bnd_Box bbox;
    BRepBndLib::Add(shape, bbox);

    if (bbox.IsVoid()) {
        return 0.0;
    }

    gp_Pnt minCorner = bbox.CornerMin();
    gp_Pnt maxCorner = bbox.CornerMax();

    return minCorner.Distance(maxCorner);
}

double GeometryUtils::ComputeBoundingBoxDiagonal(const std::vector<TopoDS_Shape>& shapes)
{
    if (shapes.empty()) {
        return 0.0;
    }

    Bnd_Box bbox;
    for (const auto& shape : shapes) {
        if (!shape.IsNull()) {
            BRepBndLib::Add(shape, bbox);
        }
    }

    if (bbox.IsVoid()) {
        return 0.0;
    }

    gp_Pnt minCorner = bbox.CornerMin();
    gp_Pnt maxCorner = bbox.CornerMax();

    return minCorner.Distance(maxCorner);
}
