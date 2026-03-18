//
// Created by Gyp Gary on 2025/1/8.
//

#ifndef TOPO_VERTEX_H
#define TOPO_VERTEX_H

#include <iostream>
#include "../geom/Point.h"

using namespace std;

class Vertex {
public:
    Vertex() = default;
    ~Vertex() = default;
    Vertex(const Point& point);
    Vertex(double x, double y, double z);

    // 重载 == 运算符，以便在 map 中使用
    bool operator==(const Vertex& other) const;
    bool operator<(const Vertex& other) const;

    double operator[](int index);  // 非常量版本
    const double operator[](int index) const;  // 常量版本，给const对象使用

    void print() const;

    Point point;  // 几何信息

//    // 自定义哈希函数
//    struct Hash {
//        size_t operator()(const Vertex& v) const {
//            size_t h1 = std::hash<int>{}(v.point[0]);
//            size_t h2 = std::hash<int>{}(v.point[1]);
//            size_t h3 = std::hash<int>{}(v.point[2]);
//            return h1 ^ (h2 << 1) ^ (h3 << 2);
//        }
//    };
};

// 定义Vertex的哈希函数
namespace std {
    template<>
    struct hash<Vertex> {
        size_t operator()(const Vertex& v) const{
            size_t h1 = std::hash<double>{}(v.point[0]);
            size_t h2 = std::hash<double>{}(v.point[1]);
            size_t h3 = std::hash<double>{}(v.point[2]);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}


#endif //TOPO_VERTEX_H