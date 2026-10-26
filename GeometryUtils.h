#pragma once
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <string>
#include <vector>

class GeometryUtils {
public:
	// 将归一化参数限制到 [0, 1]；不作为有限值有效性检查。
	static double Clamp01(double value);
	// 点的三个坐标均为有限值时返回 true。
	static bool IsFinitePoint(const gp_Pnt& point);

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

	/// <summary>
	/// 判断边在指定面上是否已有可用 pcurve。
	/// </summary>
	static bool HasValidPCurveOnFace(const TopoDS_Edge& edge, const TopoDS_Face& face);
	static bool HasValidPCurveOnFace(const TopoDS_Edge& edge,
	                                 const TopoDS_Face& face,
	                                 double tolerance,
	                                 std::string* errorMessage = nullptr);

	/// <summary>
	/// 为边在指定面上补建并校验 pcurve，同时修正 2D 方向和 same-parameter。
	/// 常规修复失败时，对非周期、非闭合支撑面尝试有界自适应投影插值重建。
	/// </summary>
	static bool RepairPCurveForEdgeOnFace(TopoDS_Edge& edge,
	                                      const TopoDS_Face& face,
	                                      double tolerance,
	                                      std::string* errorMessage = nullptr);

	/// <summary>
	/// 针对引用该边的所有面补建并校验 pcurve。
	/// </summary>
	static bool RepairPCurvesForEdge(TopoDS_Edge& edge,
	                                 const TopTools_ListOfShape& faces,
	                                 double tolerance,
	                                 std::string* errorMessage = nullptr);
};
