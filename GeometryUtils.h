#pragma once
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Vec.hxx>

class GeometryUtils {
public:
	static double ComputeEdgeLength(const TopoDS_Edge& edge);
	/// <summary>
	/// 检查这条边与相邻边的夹角，如果夹角过小，说明这条边可能是一个尖锐的边，不适合被合并
	/// </summary>
	/// <returns></returns>
	static double ComputeIntersectionAngle(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex commonV);

	
};