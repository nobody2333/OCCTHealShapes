#pragma once
#include<TopAbs_ShapeEnum.hxx>
#include <TopoDS_Vertex.hxx>
#include<TopTools_IndexedDataMapOfShapeListofShape.hxx>
#include<TopExp.hxx>

class TopologyUtils {
public:
	/// <summary>
	/// 构建顶点-边映射
	/// </summary>
	/// <param name="shape">:需要构建的拓扑形状</param>
	/// <param name="map">:构建完成的映射表</param>
	static void BuildVertexEdgeMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map);
	/// <summary>
	/// 获取顶点的度数
	/// </summary>
	/// <param name="vertex">:需要获取度的顶点</param>
	/// <param name="map">:先前构建的映射表</param>
	/// <returns></returns>
	static int GetVertexDegree(const TopoDS_Vertex& vertex, const TopTools_IndexedDataMapOfShapeListOfShape& map);
};