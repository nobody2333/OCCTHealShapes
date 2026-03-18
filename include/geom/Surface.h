#ifndef GEOM_SURFACE_H
#define GEOM_SURFACE_H


#include "TinyVector.h"
#include "Point.h"

/**
 * @brief 几何曲面基类
 */
class Surface
{
public:
    Surface() = default;
    virtual ~Surface() = default;

    virtual Vec4d getUVScale() = 0;  // return umin, umax, vmin, vmax;

    virtual Vec6d getXYZScale() = 0;  // return xmin, ymin, zmin, xmax, ymax, zmax

    virtual int get_u_control_points_num() {
        return -1;
    };

    virtual int get_v_control_points_num() {
        return -1;
    };

    virtual Vec3d get_control_point(int u_index, int v_index) {
        return { 0, 0, 0 };
    };

    // 下面四个均被换成了=0
    virtual int get_points_num() {
        return -1;
    };

    virtual int get_elements_num() {
        return -1;
    };

    virtual Point get_point(int index) {
        return { 0, 0, 0 };
    };

    virtual int get_point_index(int index) {
        return -1;
    };

    virtual bool isDiscrete() = 0;  // 判断曲面是连续的还是离散的，连续的返回false、离散的返回true

    //virtual void set_elements_num(int n) {}

    /**
     * @brief given a coord from uv space, return the cooresponding coord in 3d space
     *
     * @param [in] uv {u,v}
     * @return std::array<double, 3> {x,y,z}
     */
    virtual Vec3d eval(double u, double v) = 0;

    /**
     * @brief
     *
     * @param[in] uv uv参数平面坐标
     * @param[out] du u方向导数
     * @param[out] dv V方向导数
     */
    virtual void d1(const Vec2d& uv, Vec3d& du, Vec3d& dv) = 0;

    /**
     * @brief
     *
     * @param[in] uv uv参数平面坐标
     * @param[out] duu u方向二阶导数
     * @param[out] dvv v方向二阶导数
     * @param[out] duv uv方向二阶导数
     */
    virtual void d2(const Vec2d& uv, Vec3d& duu, Vec3d& dvv, Vec3d& duv) = 0;
	
	
    /**
     * @brief 
	 Evaluate the unitary vectors normal to the surface the input (u,v) parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     * @param[out] norm 曲线在U方向的法向。
     */
    virtual void normal(double u, double v, Vec3d& norm) = 0;
	
	
	/**
     * @brief 
	 Evaluate the mean and gaussian curvatures of the surface the input (u,v) parametrization
     *
     * @param[in] u  输入的UV参数空间中的U值。
     */
    virtual Vec2d curvature(double u, double v) = 0;


    /**
     * @brief
     *
     * @param xyz 三维空间坐标
     * @param eps 投影误差阈值
     * @param hint 提示值
     * @param res_uv 投影结果(u, v)坐标
     */
//    virtual void project(const Vec3d& xyz, double eps, const Vec2d& hint, Vec2d& resUV) = 0;


    /**
     * @brief
     *
     * @param xyz 三维空间坐标
     * @param u 投影结果 u 坐标
     * @param v 投影结果 v 坐标
     * @param prj 投影结果三维坐标
     */
    virtual void project(const Vec3d& xyz, double &u, double &v, Vec3d& prj) = 0;
};


#endif