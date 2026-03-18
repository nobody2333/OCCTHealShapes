//
// Created by Gyp Gary on 2025/1/10.
//

#ifndef GEOM_POINT_H
#define GEOM_POINT_H

#include <iostream>
#include "TinyVector.h"

using namespace std;

class Point {
public:
    Point() = default;
    ~Point() = default;

    Point(double x, double y, double z);
//    Point(Vec3d &point);
    Point(const Point &point);

    void output() const;  // 打印函数

    double& operator[](int index);  // 非常量版本
    const double& operator[](int index) const;  // 常量版本，给const对象使用的

private:
    Vec3d point;  // 成员变量，三维顶点
};

#endif //GEOM_POINT_H