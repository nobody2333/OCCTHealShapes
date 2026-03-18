//
// Created by Gyp Gary on 2025/4/1.
//
#pragma once

#include <iostream>
#include <vector>
#include "TinyVector.h"

using namespace std;

// 这是一个面片
class GeomElement {
public:
	GeomElement() = default;
	~GeomElement() = default;

	GeomElement(vector<int> points_index, int points_num);

	int get_points_num();
	void set_points_num(int n);

	int get_points_index(int i);
	void set_points_index(int n);

private:
	vector<int> _points_index;  // 面片顶点的索引下标
	int _points_num;            // 每个面片的顶点数量
};