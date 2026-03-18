//
// Created by Gyp Gary on 2025/1/8.
//
#pragma once

#include "Curve.h"
#include "GeomElement.h"
#include <iostream>

using namespace std;

class DiscreteCurve : public Curve {
public:
	// 这里是离散曲面，因此返回true
	bool isDiscrete() override;

private:
	vector<GeomElement> _elements;  // 几何单元
};