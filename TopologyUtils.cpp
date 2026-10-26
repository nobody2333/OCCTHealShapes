#include "TopologyUtils.h"

#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>

void TopologyUtils::BuildVertexEdgeMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map)
{
    map.Clear();
    if (shape.IsNull()) {
        return;
    }

    TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, map);
    RemoveDuplicateAncestors(map);
}

void TopologyUtils::BuildEdgeFaceMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map)
{
    map.Clear();
    if (shape.IsNull()) {
        return;
    }

    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, map);
    RemoveDuplicateAncestors(map);
}

void TopologyUtils::BuildFaceShellMap(
    const TopoDS_Shape& shape,
    TopTools_IndexedDataMapOfShapeListOfShape& map)
{
    map.Clear();
    if (shape.IsNull()) {
        return;
    }

    TopExp::MapShapesAndAncestors(shape, TopAbs_FACE, TopAbs_SHELL, map);
    RemoveDuplicateAncestors(map);
}

int TopologyUtils::GetVertexDegree(const TopoDS_Vertex& vertex, const TopTools_IndexedDataMapOfShapeListOfShape& map)
{
    if (map.Contains(vertex)) {
        return map.FindFromKey(vertex).Extent();
    }
    return 0;
}

int TopologyUtils::CountFreeEdges(const TopoDS_Shape& shape)
{
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    BuildEdgeFaceMap(shape, edgeFaceMap);

    int count = 0;
    for (int i = 1; i <= edgeFaceMap.Extent(); ++i) {
        if (edgeFaceMap.FindFromIndex(i).Extent() <= 1) {
            ++count;
        }
    }
    return count;
}

const char* TopologyUtils::ShapeTypeName(TopAbs_ShapeEnum shapeType)
{
    switch (shapeType) {
    case TopAbs_COMPOUND:
        return "COMPOUND";
    case TopAbs_COMPSOLID:
        return "COMPSOLID";
    case TopAbs_SOLID:
        return "SOLID";
    case TopAbs_SHELL:
        return "SHELL";
    case TopAbs_FACE:
        return "FACE";
    case TopAbs_WIRE:
        return "WIRE";
    case TopAbs_EDGE:
        return "EDGE";
    case TopAbs_VERTEX:
        return "VERTEX";
    case TopAbs_SHAPE:
        return "SHAPE";
    default:
        return "UNKNOWN";
    }
}

const char* TopologyUtils::ShapeTypeName(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) {
        return "NULL";
    }

    return ShapeTypeName(shape.ShapeType());
}

void TopologyUtils::RemoveDuplicateAncestors(TopTools_IndexedDataMapOfShapeListOfShape& map)
{
    for (int i = 1; i <= map.Extent(); ++i) {
        TopTools_ListOfShape& originalShapes = map.ChangeFromIndex(i);

        TopTools_MapOfShape uniqueShapeFilter;
        TopTools_ListOfShape uniqueShapes;

        for (TopTools_ListIteratorOfListOfShape it(originalShapes); it.More(); it.Next()) {
            const TopoDS_Shape& shape = it.Value();
            if (uniqueShapeFilter.Add(shape)) {
                uniqueShapes.Append(shape);
            }
        }

        originalShapes = uniqueShapes;
    }
}

bool TopologyUtils::GetEdgeVertices(const TopoDS_Edge& edge, TopoDS_Vertex& first, TopoDS_Vertex& last)
{
    first.Nullify();
    last.Nullify();
    if (edge.IsNull()) {
        return false;
    }

    TopExp::Vertices(edge, first, last, Standard_True);
    return !first.IsNull() && !last.IsNull();
}

TopoDS_Edge TopologyUtils::OrientEdgeFromTo(const TopoDS_Edge& edge, const TopoDS_Vertex& from, const TopoDS_Vertex& to)
{
    TopoDS_Edge oriented = edge;
    if (edge.IsNull() || from.IsNull() || to.IsNull()) {
        return oriented;
    }

    TopoDS_Vertex first;
    TopoDS_Vertex last;
    if (!TopologyUtils::GetEdgeVertices(oriented, first, last)) {
        return oriented;
    }

    if (first.IsSame(from) && last.IsSame(to)) {
        return oriented;
    }

    if (first.IsSame(to) && last.IsSame(from)) {
        oriented.Reverse();
    }
    return oriented;
}

bool TopologyUtils::IsSameFaceInVector(const std::vector<TopoDS_Face>& faces, const TopoDS_Face& face)
{
    for (const TopoDS_Face& existingFace : faces) {
        if (existingFace.IsSame(face)) {
            return true;
        }
    }
    return false;
}
