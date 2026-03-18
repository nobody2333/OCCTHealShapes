//
// Created by Gyp Gary on 2025/1/8.
//

#ifndef TOPO_LOOP_H
#define TOPO_LOOP_H


#include <iostream>
#include "Edge.h"
#include "Face.h"
#include <map>

using namespace std;

class Edge;

class Loop {
public:
    Loop() = default;
    ~Loop() = default;

    void addEdge(Edge& edge);  // 给loop增加一条边
    bool isOpen();  // 检查loop是否是封闭的、是封闭的则返回false，是open的则返回true

    void print() const;

    vector<Edge*> edges = {};  // loop包含的边
};


#endif //TOPO_LOOP_H