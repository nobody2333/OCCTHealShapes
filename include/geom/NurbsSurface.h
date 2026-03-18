//
// Created by Gyp Gary on 2024/10/23.
//

#ifndef GEOM_NURBSSURFACE_H
#define GEOM_NURBSSURFACE_H

#include "Surface.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <nlopt.hpp>
#include "Point.h"

using namespace std;

class NurbsSurface : public Surface {
public:
    NurbsSurface() = default;
    ~NurbsSurface() = default;

//    NurbsSurface(const NurbsSurface &nurbsSurface);

    NurbsSurface(const vector<vector<Vec3d>> &control_points, const vector<double> weights, const int u_order, const int v_order);

    NurbsSurface(const vector<vector<Vec3d>> &control_points, const vector<double> weights, const int u_order, const int v_order, const vector<double>& u_knots, const vector<double>& v_knots);

    // 初始化函数，用于初始化曲线的控制点、权重、阶数
    void initialize(const vector<vector<Vec3d>> &control_points, const vector<double> weights, const int u_order, const int v_order);

    // 初始化函数，用于初始化曲线的控制点、权重、阶数和节点向量
    void initialize(const vector<vector<Vec3d>> &control_points, const vector<double> weights, const int u_order, const int v_order, const vector<double>& u_knots, const vector<double>& v_knots);

    // 1.13
    vector<vector<Vec3d>> get_control_points() const {
        return this -> control_points;
    }
    vector<double> get_weights() const {
        return this -> weights;
    }
    int get_u_order() const {
        return this -> u_order;
    }
    int get_v_order() const {
        return this -> v_order;
    }
    vector<double> get_u_knots() const {
        return this -> u_knots;
    }
    vector<double> get_v_knots() const {
        return this -> v_knots;
    }


    vector<double> generate_knots(int n, int m);

    // 生成开放式节点向量
    void knot(int n, vector<double> knots, int m);

    // 生成周期性结点向量
    void knotu(int n, vector<double> knots, int m);

    Vec4d getUVScale() override; // return umin, umax, vmin, vmax;

    Vec6d getXYZScale() override; // return xmin, ymin, zmin, xmax, ymax, zmax

    int get_u_control_points_num() override;
    int get_v_control_points_num() override;

    Vec3d get_control_point(int u_index, int v_index) override;

    //int get_points_num();

    //int get_elements_num();

    //Point get_point(int index);

    //int get_point_index(int index);

    // 判断该曲面是否是离散的，本类为连续的NurbsSurface，是离散的，因此返回false
    bool isDiscrete() override;

    /**
     * @brief 计算基函数值
     *
     * @param[in] t 输入的参数值
     * @param[in] n u或v维度的控制顶点个数
     * @param[in] m 基函数的阶数
     * @param[in] knots 节点向量
     */
    vector<double> basic_func(double t, int n, int m, vector<double> knots);

    // 计算所有的基函数与权重的和，作为分母
    double sum_basis_weights(const vector<double> &u_basis, const vector<double> &v_basis);

    /**
     * @brief given a coord from uv space, return the cooresponding coord in 3d space
     *
     * @param [in] uv {u,v}
     * @return std::array<double, 3> {x,y,z}
     */
    Vec3d eval(double u, double v) override;

    // 计算基函数关于参数t的一阶导数
    // t:参数、p:阶数、n：该维度控制点的数量、knots:节点向量
    vector<double> basic_func_deriv(double t, int p, int n, vector<double> knots);

    /**
     * @brief
     *
     * @param[in] uv uv参数平面坐标
     * @param[out] du u方向导数
     * @param[out] dv V方向导数
     */
    void d1(const Vec2d& uv, Vec3d& du, Vec3d& dv) override;

    // 计算基函数关于参数t的一阶导数
    // t:参数、p:阶数、n：该维度控制点的数量、knots:节点向量
    vector<double> basic_func_second_deriv(double t, int p, int n, vector<double> knots);

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
     * @param[in] v  输入的UV参数空间中的V值。
     * @param[out] norm 计算给定参数 u 和 v 的NURBS曲面上的单位法向量
     */
    void normal(double u, double v, Vec3d& norm) override;

    /**
     * @brief
     Evaluate the mean and gaussian curvatures of the surface the input (u,v) parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[in] v  输入的UV参数空间中的V值。
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
    void project(const Vec3d& xyz, double &u, double &v, Vec3d& prj) override;

    // 目标函数
    static double objectiveFunctionWrapper(const std::vector<double>& params, std::vector<double>& grad, void* data);

private:
    vector<double> u_knots;
    vector<double> v_knots;

    // 类似“二维数组”，[i][j]中i是u维度的，j是v维度的
    vector<vector<Vec3d>> control_points;
    vector<double> weights;

    int u_order;  // u方向曲线的阶数
    int v_order;  // v方向曲线的阶数
};


#endif //GEOM_NURBSSURFACE_H