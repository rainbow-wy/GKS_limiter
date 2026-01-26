/**
 * @file riemann_problem_subcell_1d.h
 * @brief 1D子单元GKS求解器测试 - 黎曼问题
 * 
 * 使用方法：在 riemann_problem_subcell_1d.cpp 中选择想要运行的算例
 */

#pragma once
#include "gks_subcell_1d_solver.h"

/**
 * @brief 使用子单元GKS求解器计算1D黎曼问题
 * 
 * 算例配置在 cpp 文件中定义，通过注释切换
 */
void riemann_problem_subcell_1d();

/**
 * @brief 1D黎曼问题初始条件设置
 * 
 * @param fluids      [out] 流场数组
 * @param rho_L       [in]  左侧密度
 * @param u_L         [in]  左侧速度
 * @param p_L         [in]  左侧压力
 * @param rho_R       [in]  右侧密度
 * @param u_R         [in]  右侧速度
 * @param p_R         [in]  右侧压力
 * @param x_discon    [in]  间断位置
 * @param block       [in]  网格块信息
 */
void ICfor1dSubcell(
    Fluid1d* fluids,
    double rho_L, double u_L, double p_L,
    double rho_R, double u_R, double p_R,
    double x_discon,
    Block1d& block
);
