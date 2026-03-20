#include "GeometryUtils.h"

#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Extrema_ExtPC.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>

#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double GeometryUtils::ComputeEdgeLength(const TopoDS_Edge& edge)
{
    if (edge.IsNull()) {
        return 0.0;
    }

    // 获取边的几何曲线
    double firstParam, lastParam;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, firstParam, lastParam);

    if (curve.IsNull()) {
        // 如果没有曲线（退化边），计算两个顶点之间的直线距离
        TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
        TopoDS_Vertex v2 = TopExp::LastVertex(edge);
        
        if (v1.IsNull() || v2.IsNull()) {
            return 0.0;
        }

        gp_Pnt p1 = BRep_Tool::Pnt(v1);
        gp_Pnt p2 = BRep_Tool::Pnt(v2);
        
        return p1.Distance(p2);
    }

    // 创建曲线适配器并计算长度
    GeomAdaptor_Curve curveAdaptor(curve, firstParam, lastParam);
    return GCPnts_AbscissaPoint::Length(curveAdaptor, firstParam, lastParam);
}

double GeometryUtils::ComputeIntersectionAngle(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex commonV)
{
    if (e1.IsNull() || e2.IsNull() || commonV.IsNull()) {
        return -1.0; // 返回负值表示错误
    }

    // 获取两条边的几何曲线
    TopLoc_Location loc1, loc2;
    double first1, last1, first2, last2;
    
    Handle(Geom_Curve) curve1 = BRep_Tool::Curve(e1, loc1, first1, last1);
    Handle(Geom_Curve) curve2 = BRep_Tool::Curve(e2, loc2, first2, last2);

    if (curve1.IsNull() || curve2.IsNull()) {
        return -1.0; // 无法获取曲线
    }

    // 应用位置变换
    if (!loc1.IsIdentity()) {
        curve1 = Handle(Geom_Curve)::DownCast(curve1->Transformed(loc1.Transformation()));
    }
    if (!loc2.IsIdentity()) {
        curve2 = Handle(Geom_Curve)::DownCast(curve2->Transformed(loc2.Transformation()));
    }

    // 获取公共顶点的坐标（BRep_Tool::Pnt()直接返回全局坐标）
    gp_Pnt commonPoint = BRep_Tool::Pnt(commonV);

    // 计算在公共顶点处的参数值
    double param1, param2;
    
    // 确定 e1 在公共顶点处的参数
    TopoDS_Vertex v1First = TopExp::FirstVertex(e1);
    TopoDS_Vertex v1Last = TopExp::LastVertex(e1);
    
    if (commonV.IsSame(v1First)) {
        param1 = first1;
    } else if (commonV.IsSame(v1Last)) {
        param1 = last1;
    } else {
        // 公共顶点不是 e1 的端点，尝试投影
        GeomAdaptor_Curve adaptor1(curve1, first1, last1);
        Extrema_ExtPC proj1(commonPoint, adaptor1);
        if (proj1.IsDone() && proj1.NbExt() > 0) {
            param1 = proj1.Point(1).Parameter();
        } else {
            return -1.0;
        }
    }

    // 确定 e2 在公共顶点处的参数
    TopoDS_Vertex v2First = TopExp::FirstVertex(e2);
    TopoDS_Vertex v2Last = TopExp::LastVertex(e2);
    
    if (commonV.IsSame(v2First)) {
        param2 = first2;
    } else if (commonV.IsSame(v2Last)) {
        param2 = last2;
    } else {
        // 公共顶点不是 e2 的端点，尝试投影
        GeomAdaptor_Curve adaptor2(curve2, first2, last2);
        Extrema_ExtPC proj2(commonPoint, adaptor2);
        if (proj2.IsDone() && proj2.NbExt() > 0) {
            param2 = proj2.Point(1).Parameter();
        } else {
            return -1.0;
        }
    }

    // 计算两条曲线在公共顶点处的切向量
    gp_Vec tangent1, tangent2;
    
    try {
		// 计算一阶导数（切向量）,方向是参数值从小到大增长的方向
        gp_Pnt pointOnCurve1;
        curve1->D1(param1, pointOnCurve1, tangent1);
        
        gp_Pnt pointOnCurve2;
        curve2->D1(param2, pointOnCurve2, tangent2);
    } catch (...) {
        return -1.0;
    }

    // 确保切向量从公共顶点向外指
    // 如果切向量指向参数减小的方向，则反向
    if (commonV.IsSame(v1First)) {
        // 在起点处，导数方向是从起点到终点，已经是向外的
        // 但如果参数是反向的（first > last），需要反向
        if (first1 > last1) {
            tangent1.Reverse();
        }
    } else if (commonV.IsSame(v1Last)) {
        // 在终点处，导数方向是从起点到终点，需要反向
        tangent1.Reverse();
        if (first1 > last1) {
            tangent1.Reverse(); // 再反向回来
        }
    }

    if (commonV.IsSame(v2First)) {
        if (first2 > last2) {
            tangent2.Reverse();
        }
    } else if (commonV.IsSame(v2Last)) {
        tangent2.Reverse();
        if (first2 > last2) {
            tangent2.Reverse();
        }
    }

    // 标准化切向量
    double mag1 = tangent1.Magnitude();
    double mag2 = tangent2.Magnitude();

    if (mag1 < 1e-10 || mag2 < 1e-10) {
        return -1.0; // 切向量太小，无法计算角度
    }

    tangent1.Normalize();
    tangent2.Normalize();

    // 计算点积并求夹角
    double dotProduct = tangent1.Dot(tangent2);
    
    // 限制点积范围在 [-1, 1] 之间，避免浮点误差
    if (dotProduct > 1.0) dotProduct = 1.0;
    if (dotProduct < -1.0) dotProduct = -1.0;

    // 计算夹角（弧度）
    double angle = std::acos(dotProduct);

    return angle;
}
