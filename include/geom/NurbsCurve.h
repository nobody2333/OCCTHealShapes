//
// Created by Gyp Gary on 2024/10/23.
//

#ifndef GEOM_NURBSCURVE_H
#define GEOM_NURBSCURVE_H

#include "Curve.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <nlopt.hpp>

using namespace std;

class NurbsCurve : public Curve {
public:
    NurbsCurve() = default;
    ~NurbsCurve() = default;

    // 构造函数
    NurbsCurve(const vector<Vec3d>& control_points, const vector<double> weights, const int order);

    NurbsCurve(const vector<Vec3d>& control_points, const vector<double> weights, const int order, const vector<double>& knots);

    // 初始化函数，用于初始化曲线的控制点、权重、阶数
    void initialize(const vector<Vec3d>& control_points, const vector<double> weights, const int order);

    // 初始化函数，用于初始化曲线的控制点、权重、阶数和节点向量
    void initialize(const vector<Vec3d>& control_points, const vector<double> weights, const int order, const vector<double>& knots);

    // 生成节点向量
    vector<double> generate_knots();

    // 生成开放式节点向量
    void knot();
    // 生成周期性结点向量
    void knotu();

    // 返回控制点数量
    int get_control_points_num() override;

    // 根据指定索引获取控制点
    Vec3d get_control_point(int index) override;

    // 判断是离散的还是连续的，这里是NurbsCurve连续曲线，是连续的，因此返回false
    bool isDiscrete() override;

    /**
     * @brief 计算基函数值
     *
     * @param[in] u 输入的参数值
     */
    vector<double> basic_func(double u);

    /**
     * @brief 计算基函数值
     *
     * @param[in] u 输入的参数值
     * @param[in] k 基函数的阶数
     */
    vector<double> basic_func(double u, int m);

    /**
     * @brief 用于获取在指定U值处，曲线在U点的坐标值。
     *
     * @param [in] u 输入的UV参数空间中的U值。
     * @return Vec3 曲线在U点的坐标值.
     */
    Vec3d eval(double u) override;


    // 辅助函数：计算 NURBS 基函数的一阶导数
    vector<double> basic_func_deriv(double u);

    vector<double> basic_func_deriv(double u, int p);

    /**
     * @brief
     *
     * @param[in] uv 输入的UV参数空间中的U值。
     * @param[out] du u方向导数
     */
    void d1(double u, Vec3d& du) override;

    // 辅助函数：计算 NURBS 基函数的二阶导数
    vector<double> basic_func_second_deriv(double u);

    vector<double> basic_func_second_deriv(double u, int p);

    /**
     * @brief 获取在指定U值处，曲线在U方向的二阶导数。
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] d2u 曲线在U方向的二阶导数。
     */
    void d2(double u, Vec3d& d2u) override;

    /**
     * @brief
	 Evaluate the unitary normal vector to the curve for the input u-parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] norm 曲线在U方向的法向。
     */
    void normal(double u, Vec3d& norm) override;

    /**
     * @brief
     Evaluate the curvature of the curve for the input u-parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     */
    double curvature(double u) override;

    /**
     * @brief 获取曲线总长
     *
     * @return[out] 曲线总长
     */
    double curve_length() override;

    /**
     * @brief 获取起点和终点分别为参数坐标u1,u2间曲线的弧长，并返回
     *
     * @return[out] arc(u1, u2)弧长
     */
    double curve_arc_length(double u1, double u2) override;
    double curve_arc_length1(double u1, double u2);
    double curve_arc_length2(double u1, double u2);
    double curve_arc_length3(double u1, double u2);

    // 给定高斯勒让德节点和权重数量n，在给定n的情况下计算
    double curve_arc_length(double u1, double u2, int n);

    // 辅助函数：计算弧长函数的导数
    double arc_length_deriv(double u);

    /**
     * @brief 获取从参数坐标u1开始弧长长度为length的终点参数坐标，并返回
     *
     * @return[out] 终点参数坐标
     */
    double curve_arc_length_u(double length, double u1) override;

    /**
     * @brief 将3D空间中的一个点投影到曲线上，并返回其对应的U值
     *
     * @param[in] xyz 三维空间坐标
     * @param[in] eps 投影误差阈值
     * @param[in] hint 提供的初值提示，有助于提高投影的计算效率和成功率。
     * @param[out] res_uv 投影结果(u, v)坐标
     */
    // 投影：给定一个三维空间中的点xyz，返回参数平面上的u值，使得u值对应的点prj与xyz的距离最小
    void project(const Vec3d& xyz, double& u, Vec3d& prj) override;

    /**
     * @brief 将3D空间中的一个点投影到曲线上，并返回其对应的U值
     *
     * @param[in] xyz 三维空间坐标
     * @param[in] eps 投影误差阈值
     * @param[in] hint 提供的初值提示，有助于提高投影的计算效率和成功率。
     * @param[out] res_uv 投影结果(u, v)坐标
     */
//    void project(const Vec3d& xyz, double eps, double hint, double& u) override;

    static double objectiveFunctionWrapper(const std::vector<double>& params, std::vector<double>& grad, void* data);

private:
    vector<Vec3d> control_points;
    vector<double> weights;
    int order;
    vector<double> knots;
};

#endif //GEOM_NURBSCURVE_H