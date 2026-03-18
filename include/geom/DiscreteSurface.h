//
// Created by Gyp Gary on 2025/1/8.
//
#pragma once

#include "Surface.h"
#include "GeomElement.h"
#include "Point.h"
#include <vector>
#include <iostream>

using namespace std;

class DiscreteSurface : public Surface {
public:
	DiscreteSurface() = default;
	~DiscreteSurface() = default;

    DiscreteSurface(int& elements_num, int& points_num, vector<Point>& points, vector<int>& points_index, vector<int>& element_points_num);

	//DiscreteSurface(vector<GeomElement> elements, int element_num);

    Vec4d getUVScale() override;  // return umin, umax, vmin, vmax;

    Vec6d getXYZScale() override;  // return xmin, ymin, zmin, xmax, ymax, zmax

    //int get_u_control_points_num() override;
    //int get_v_control_points_num() override;

    //Vec3d get_control_point(int u_index, int v_index) override;

    int get_points_num() override;

    int get_elements_num() override;

    Point get_point(int index) override;

    int get_point_index(int index) override;

    // 当前曲面是离散的，因此返回true
    bool isDiscrete() override;

    //void set_elements_num(int n) override;

    /**
     * @brief given a coord from uv space, return the cooresponding coord in 3d space
     *
     * @param [in] uv {u,v}
     * @return std::array<double, 3> {x,y,z}
     */
    Vec3d eval(double u, double v) override;

    /**
     * @brief
     *
     * @param[in] uv uv参数平面坐标
     * @param[out] du u方向导数
     * @param[out] dv V方向导数
     */
    void d1(const Vec2d& uv, Vec3d& du, Vec3d& dv) override;

    /**
     * @brief
     *
     * @param[in] uv uv参数平面坐标
     * @param[out] duu u方向二阶导数
     * @param[out] dvv v方向二阶导数
     * @param[out] duv uv方向二阶导数
     */
    void d2(const Vec2d& uv, Vec3d& duu, Vec3d& dvv, Vec3d& duv) override;


    /**
     * @brief
     Evaluate the unitary vectors normal to the surface the input (u,v) parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] norm 曲线在U方向的法向。
     */
    void normal(double u, double v, Vec3d& norm) override;


    /**
     * @brief
     Evaluate the mean and gaussian curvatures of the surface the input (u,v) parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     */
    Vec2d curvature(double u, double v) override;

    /**
    * @brief
    *
    * @param xyz 三维空间坐标
    * @param u 投影结果 u 坐标
    * @param v 投影结果 v 坐标
    * @param prj 投影结果三维坐标
    */
    void project(const Vec3d& xyz, double& u, double& v, Vec3d& prj) override;
     
private:
	vector<GeomElement> _elements;  // 几何面片单元
	int _elements_num;			   // 几何面片单元的数量
    int _points_num;                // 几何顶点总数
    vector<Point> _points;          // 顶点坐标数组
    vector<int> _points_index;      // 面顶点索引，描述多边形顶点连接顺序
    vector<int> _element_points_num;// 每个面的顶点数(3-三角形、4-四边形，按照目前策略是三角形)
};