//
// Created by Gyp Gary on 2025/1/8.
//

#ifndef TOPO_BODY_H
#define TOPO_BODY_H

#include <iostream>
#include <map>
#include "Face.h"

class Edge;

class Body {
public:
    Body() = default;
    ~Body() = default;

    void addFace(Face& face);
    bool is_open();  // 判断是否开放（封闭），封闭的返回false，开放的返回true

    void print() const;

public:
    vector<Face*> faces = {};  // 体包含的面
    bool isSolid = true;
};


#endif //TOPO_BODY_H