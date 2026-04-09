#include "TopologyUtils.h"

void TopologyUtils::BuildVertexEdgeMap(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeListOfShape& map) {
    map.Clear();

    TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, map);

    for (int i = 1; i <= map.Extent(); ++i) {
        //使用 ChangeFromIndex 获取原始列表的允许修改的引用
        TopTools_ListOfShape& originalEdges = map.ChangeFromIndex(i);

        TopTools_MapOfShape uniqueEdgeFilter;  //利用哈希进行快速去重
        TopTools_ListOfShape uniqueEdges;      //存放真正不重复的边

        //遍历当前顶点原本的边列表
        for (TopTools_ListIteratorOfListOfShape it(originalEdges); it.More(); it.Next()) {
            const TopoDS_Shape& edge = it.Value();

            //如果是第一次遇到这个形状，会将其加入哈希表并返回true
            if (uniqueEdgeFilter.Add(edge)) {
                uniqueEdges.Append(edge);
            }
        }

        originalEdges = uniqueEdges;
    }
}

int TopologyUtils::GetVertexDegree(const TopoDS_Vertex& vertex, const TopTools_IndexedDataMapOfShapeListOfShape& map) {
	if(map.Contains(vertex))
		return map.FindFromKey(vertex).Extent();
	return 0;
}