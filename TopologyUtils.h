#pragma once

#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <vector>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

class TopologyUtils {
public:
    // 获取考虑 edge orientation 的起止顶点；先清空输出，缺少任一端点时返回 false。
    static bool GetEdgeVertices(
        const TopoDS_Edge& edge, TopoDS_Vertex& first, TopoDS_Vertex& last);

    // 按拓扑端点调整 edge use 方向，不修改底层 TShape。
    // 输入为空、缺少端点或端点不匹配时原样返回，不表示定向成功。
    static TopoDS_Edge OrientEdgeFromTo(
        const TopoDS_Edge& edge, const TopoDS_Vertex& from, const TopoDS_Vertex& to);

    // 使用 IsSame 比较 TShape 和 location，忽略 orientation，不按几何位置合并。
    static bool IsSameFaceInVector(
        const std::vector<TopoDS_Face>& faces, const TopoDS_Face& face);

    /// <summary>
    /// 构建 vertex -> adjacent edges 映射，并去除重复 edge。
    /// </summary>
    /// <param name="shape">需要构建映射的拓扑形状。</param>
    /// <param name="map">构建完成的映射表。</param>
    static void BuildVertexEdgeMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map);

    /// <summary>
    /// 构建 edge -> adjacent faces 映射，并去除重复 face。
    /// </summary>
    /// <param name="shape">需要构建映射的拓扑形状。</param>
    /// <param name="map">构建完成的映射表。</param>
    static void BuildEdgeFaceMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map);

    /// <summary>
    /// 构建 face -> ancestor shells 映射，并去除重复 shell。
    /// </summary>
    static void BuildFaceShellMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map);

    /// <summary>
    /// 获取顶点连接的唯一 edge 数量。
    /// </summary>
    /// <param name="vertex">需要获取度数的顶点。</param>
    /// <param name="map">先前构建的 vertex -> edge 映射表。</param>
    /// <returns>顶点度数。</returns>
    static int GetVertexDegree(const TopoDS_Vertex& vertex, const TopTools_IndexedDataMapOfShapeListOfShape& map);

    /// <summary>
    /// 统计只被 0 或 1 个 face 使用的 free edge 数量。
    /// </summary>
    /// <param name="shape">需要统计的拓扑形状。</param>
    /// <returns>free edge 数量。</returns>
    static int CountFreeEdges(const TopoDS_Shape& shape);

    /// <summary>
    /// 将 OCCT shape type 枚举转换为便于日志输出的文字。
    /// </summary>
    static const char* ShapeTypeName(TopAbs_ShapeEnum shapeType);

    /// <summary>
    /// 获取 shape 的类型文字；空 shape 返回 NULL。
    /// </summary>
    static const char* ShapeTypeName(const TopoDS_Shape& shape);

private:
    /// <summary>
    /// 对 IndexedDataMap 中每个 key 对应的 shape list 做拓扑去重。
    /// </summary>
    static void RemoveDuplicateAncestors(TopTools_IndexedDataMapOfShapeListOfShape& map);
};
