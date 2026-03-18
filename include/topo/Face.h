//
// Created by Gyp Gary on 2025/1/8.
//

#ifndef TOPO_FACE_H
#define TOPO_FACE_H

#include <iostream>
#include <vector>
#include "Edge.h"
#include "Loop.h"
#include "../geom/NurbsSurface.h"
#include "../geom/Surface.h"

using namespace std;

class Edge;
class Loop;

class Face {
public:
    Face();
    ~Face() = default;

    void addEdge(Edge& edge);
    void addLoop(Loop& loop);

    void print() const;

    // 该拓扑Surface对应的曲面信息，包含连续的NurbsSurface和离散的DiscreteSurface
    Surface* geom_surface{};      // 连续曲面
    Surface* discrete_surface{};  // 离散曲面
    vector<Edge*> edges = {};     // 面包含的边
    vector<Loop*> loops = {};     // 面包含的环
    int index;                    // 面的索引
};


#endif //TOPO_FACE_H