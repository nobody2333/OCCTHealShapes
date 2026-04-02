#pragma once
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Vec.hxx>
#include <vector>

class GeometryUtils {
public:
	static double ComputeEdgeLength(const TopoDS_Edge& edge);
	/// <summary>
	/// 计算两条边在公共顶点处的切线夹角（弧度）
	/// </summary>
	/// <returns></returns>
	static double ComputeIntersectionAngle(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex commonV);

	/// <summary>
	/// 计算单个形状的包围盒对角线长度
	/// </summary>
	/// <param name="shape">输入形状</param>
	/// <returns>包围盒对角线长度，失败返回0</returns>
	static double ComputeBoundingBoxDiagonal(const TopoDS_Shape& shape);

	/// <summary>
	/// 计算多个形状的包围盒对角线长度（并集包围盒）
	/// </summary>
	/// <param name="shapes">输入形状列表</param>
	/// <returns>包围盒对角线长度，失败返回0</returns>
	static double ComputeBoundingBoxDiagonal(const std::vector<TopoDS_Shape>& shapes);

	
};