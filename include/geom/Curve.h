#ifndef GEOM_CURVE_H
#define GEOM_CURVE_H


#include "TinyVector.h"

// 定义一个曲线类型
enum CurveType
{
    Nurbs = 0,
    BSpline = 1,
    Bezier = 2
};

/**
 * @brief 几何曲线基类
 */
class Curve
{
public:
    Curve() = default;
    virtual ~Curve() = default;  // 使用默认的虚析构函数

    // 获取控制点数量
    virtual int get_control_points_num() = 0;

    // 根据指定索引获取控制点
    virtual Vec3d get_control_point(int index) = 0;

    // 判断是离散的还是连续的
    virtual bool isDiscrete() = 0;

    /**
     * @brief 用于获取在指定U值处，曲线在U点的坐标值。
     *
     * @param [in] u 输入的UV参数空间中的U值。
     * @return Vec3d 曲线在U点的坐标值.
     */
    virtual Vec3d eval(double u) = 0;

    /**
     * @brief
     *
     * @param[in] uv 输入的UV参数空间中的U值。
     * @param[out] du u方向导数
     */
    virtual void d1(double u, Vec3d& du) = 0;

    /**
     * @brief 获取在指定U值处，曲线在U方向的二阶导数。
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] d2u 曲线在U方向的二阶导数。
     */
    virtual void d2(double u, Vec3d& d2u) = 0;
	
	
    /**
     * @brief 
	 Evaluate the unitary normal vector to the curve for the input u-parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] norm 曲线在U方向的法向。
     */
    virtual void normal(double u, Vec3d& norm) = 0;
	
	/**
     * @brief 
	 Evaluate the curvature of the curve for the input u-parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     */
    virtual double curvature(double u) = 0;

    /**
     * @brief 获取曲线总长
     *
     * @return[out] 曲线总长
     */
    virtual double curve_length() = 0;
	
	/**
     * @brief 获取从参数坐标u1开始弧长长度为length的终点参数坐标，并返回
     *
     * @return[out] 终点参数坐标
     */
	virtual double curve_arc_length_u(double length, double u1) = 0;

	/**
     * @brief 获取起点和终点分别为参数坐标u1,u2间曲线的弧长，并返回
     *
     * @return[out] arc(u1, u2)弧长
     */
	virtual double curve_arc_length(double u1, double u2) = 0;
	
    /**
     * @brief 将3D空间中的一个点投影到曲线上，并返回其对应的U值
     *
     * @param[in] xyz 三维空间坐标
     * @param[in] eps 投影误差阈值
     * @param[in] hint 提供的初值提示，有助于提高投影的计算效率和成功率。
     * @param[out] res_uv 投影结果(u, v)坐标
     */
//    virtual void project(const Vec3d& xyz, double eps, double hint, double& u) = 0;

	/**
     * @brief 将3D空间中的一个点投影到曲线上，并返回其对应的U值
     *
     * @param[in] xyz 三维空间坐标
     * @param[in] eps 投影误差阈值
     * @param[in] hint 提供的初值提示，有助于提高投影的计算效率和成功率。
     * @param[out] res_uv 投影结果(u, v)坐标
     */
    virtual void project(const Vec3d& xyz, double& u, Vec3d& prj) = 0;

private:
    CurveType curveType;
};


#endif  // GEOM_CURVE_H