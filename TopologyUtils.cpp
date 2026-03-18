#include "TopologyUtils.h"

void TopologyUtils::BuildVertexEdgeMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map) {
	map.Clear();
	TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, map);
}

int TopologyUtils::GetVertexDegree(const TopoDS_Vertex& vertex, const TopTools_IndexedDataMapOfShapeListOfShape& map) {
	if(map.Contains(vertex))
		return map.FindFromKey(vertex).Extent();
	return 0;
}