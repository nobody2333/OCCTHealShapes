//
// Created by Gyp Gary on 2025/1/8.
//

#ifndef TOPO_EDGE_H
#define TOPO_EDGE_H

#include <iostream>
#include "Vertex.h"
#include "Face.h"
#include "Loop.h"
#include "../geom/NurbsCurve.h"
#include "../geom/Curve.h"

using namespace std;

class Face;
class Loop;

class Edge {
public:
    Edge() = default;
    ~Edge() = default;
    Edge(const Edge& other);
    Edge(const Vertex& v1, const Vertex& v2);
    Edge(Curve* curve);
    Edge(const Vertex& v1, const Vertex& v2, Curve* curve);

    bool operator==(const Edge& other) const;
    bool operator<(const Edge& other) const;

    void print() const;

    // 组成边的两个顶点
    Vertex vertex1;  // 起始点
    Vertex vertex2;  // 终点
    Curve* curve;    // 该拓扑Edge对应的曲线信息，包含连续NurbsCurve和离散DiscreteCurve两种
    vector<Face*> faces = {};  // 连接的面
    vector<Loop*> loops = {};  // 连接的环
};


#endif //TOPO_EDGE_H