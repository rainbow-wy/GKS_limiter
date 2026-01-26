/**
 * @file interface_1d_limiter.h
 * @brief 1D界面通量计算与保正限制
 * 
 * 宏观界面：混合高阶GKS与低阶KFVS + 试探-检查-修正
 * 子单元内部：二阶KFVS + 线性重构
 */

#pragma once
#include "subcell_1d_data.h"
#include "subcell_1d_tools.h"

//=============================================================================
// 宏观界面通量计算
//=============================================================================

/**
 * @brief 计算宏观界面的高阶GKS通量
 * 
 * 流程:
 * 1. 保存 gks1dsolver
 * 2. 设置 gks1dsolver = gks2nd
 * 3. 调用 WENO5_AO 重构
 * 4. 调用 GKS 计算通量
 * 5. 恢复 gks1dsolver
 * 
 * @param flux       [out] 高阶通量 [3]
 * @param fluids     [in]  宏观单元数组
 * @param iface_idx  [in]  界面索引 (界面在 cell[i] 和 cell[i+1] 之间, idx=i+1)
 * @param dt         [in]  时间步长
 */
void ComputeMacroFluxHigh(
    double* flux,
    Fluid1d* fluids,
    int iface_idx,
    double dt
);

/**
 * @brief 计算宏观界面的低阶KFVS通量 (一阶保正)
 * 
 * 流程:
 * 1. 保存 gks1dsolver
 * 2. 设置 gks1dsolver = kfvs1st
 * 3. 一阶重构 (直接使用单元平均值)
 * 4. 调用 GKS 计算通量
 * 5. 恢复 gks1dsolver
 * 
 * @param flux       [out] 低阶通量 [3]
 * @param fluids     [in]  宏观单元数组
 * @param iface_idx  [in]  界面索引
 * @param dt         [in]  时间步长
 */
void ComputeMacroFluxLow(
    double* flux,
    Fluid1d* fluids,
    int iface_idx,
    double dt
);

/**
 * @brief 计算混合通量
 * 
 * F_blend = (1 - α_bar) * F_high + α_bar * F_low
 * α_bar = 0.5 * (α_left + α_right)
 * 
 * @param flux_blend [out] 混合通量 [3]
 * @param flux_high  [in]  高阶通量 [3]
 * @param flux_low   [in]  低阶通量 [3]
 * @param alpha_L    [in]  左单元α
 * @param alpha_R    [in]  右单元α
 */
void BlendFlux(
    double* flux_blend,
    const double* flux_high,
    const double* flux_low,
    double alpha_L,
    double alpha_R
);

/**
 * @brief 试探性更新界面两侧子单元
 * 
 * @param u_L_new [out] 左子单元更新后状态 [3]
 * @param u_R_new [out] 右子单元更新后状态 [3]
 * @param sc_L    [in]  左宏观单元最右子单元
 * @param sc_R    [in]  右宏观单元最左子单元
 * @param flux    [in]  界面通量 [3]
 * @param dt      [in]  时间步长
 */
void TentativeUpdate(
    double* u_L_new,
    double* u_R_new,
    const SubCell1d& sc_L,
    const SubCell1d& sc_R,
    const double* flux,
    double dt
);

/**
 * @brief 计算修正系数θ (二分法)
 * 
 * F_final = θ * F_cand + (1 - θ) * F_low
 * 找最大的θ使更新后状态满足保正性
 * 
 * @param sc_L   [in] 左子单元
 * @param sc_R   [in] 右子单元
 * @param F_cand [in] 候选混合通量 [3]
 * @param F_low  [in] 低阶通量 [3]
 * @param dt     [in] 时间步长
 * @return θ ∈ [0, 1]
 */
double ComputeCorrectionTheta(
    const SubCell1d& sc_L,
    const SubCell1d& sc_R,
    const double* F_cand,
    const double* F_low,
    double dt
);

/**
 * @brief 计算宏观界面最终通量 (完整的试探-检查-修正流程)
 * 
 * @param flux_final [out] 最终通量 [3]
 * @param flux_high  [in]  高阶通量 [3]
 * @param flux_low   [in]  低阶通量 [3]
 * @param sc_L       [in]  左宏观单元最右子单元
 * @param sc_R       [in]  右宏观单元最左子单元
 * @param alpha_L    [in]  左宏观单元α
 * @param alpha_R    [in]  右宏观单元α
 * @param dt         [in]  时间步长
 */
void ComputeMacroFluxFinal(
    double* flux_final,
    const double* flux_high,
    const double* flux_low,
    const SubCell1d& sc_L,
    const SubCell1d& sc_R,
    double alpha_L,
    double alpha_R,
    double dt
);

//=============================================================================
// 子单元内部通量计算
//=============================================================================

/**
 * @brief 子单元线性重构 (VanLeer限制器)
 * 
 * 计算子单元界面处的左右值和斜率
 * 
 * @param left_val   [out] 界面左值 [3]
 * @param right_val  [out] 界面右值 [3]
 * @param left_der   [out] 左侧斜率 [3]
 * @param right_der  [out] 右侧斜率 [3]
 * @param sc_left    [in]  左子单元
 * @param sc_center  [in]  中心子单元 (用于计算斜率)
 * @param sc_right   [in]  右子单元
 */
void SubcellLinearRecon(
    double* limited_slope,
    const SubCell1d& sc_left,
    const SubCell1d& sc_center,
    const SubCell1d& sc_right
);

/**
 * @brief 计算子单元内部界面的二阶KFVS通量
 * 
 * 流程:
 * 1. 保存 gks1dsolver
 * 2. 设置 gks1dsolver = kfvs2nd
 * 3. 使用线性重构值填充 Interface1d
 * 4. 调用 GKS 计算通量
 * 5. 恢复 gks1dsolver
 * 
 * @param flux      [out] 通量 [3]
 * @param sc_left   [in]  界面左侧子单元
 * @param sc_right  [in]  界面右侧子单元
 * @param left_der  [in]  左侧斜率 [3] (来自线性重构)
 * @param right_der [in]  右侧斜率 [3]
 * @param dt        [in]  时间步长
 */
void ComputeInternalFlux(
    double* flux,
    const SubCell1d& sc_left,
    const SubCell1d& sc_right,
    const double* left_der,
    const double* right_der,
    double dt
);

//=============================================================================
// 全场通量计算
//=============================================================================

/**
 * @brief 计算所有宏观界面通量
 * 
 * @param fluids_sc [in/out] 子单元数组 (flux_left, flux_right 被更新)
 * @param fluids    [in]     宏观单元数组
 * @param block     [in]     网格块信息
 * @param dt        [in]     时间步长
 */
void ComputeAllMacroFluxes(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block,
    double dt
);

/**
 * @brief 计算所有子单元内部通量
 * 
 * @param fluids_sc [in/out] 子单元数组 (flux_internal 被更新)
 * @param block     [in]     网格块信息
 * @param dt        [in]     时间步长
 */
void ComputeAllInternalFluxes(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double dt
);

//=============================================================================
// 缩放限制器 (Scaling Limiter) - Zhang & Shu, 2010 算法
//=============================================================================

/**
 * @brief 对单个宏观单元内的子单元应用缩放限制器
 * 
 * 此函数是算法最后一道保险，在通量更新后调用。
 * 前置条件：宏观单元平均值已经保正（密度>0, 压力>0）。
 * 
 * 算法流程（串行修正策略）：
 *   1. 计算宏观平均值 U_mean（作为安全锚点）
 *   2. 第一步：密度限制 - 线性缩放使所有子单元密度 >= epsilon
 *   3. 第二步：压力限制 - 二次方程求解使所有子单元压力 >= epsilon
 * 
 * @param cell    [in/out] 宏观单元（含5个子单元），子单元值会被就地修改
 * @param epsilon [in]     安全下限阈值（默认 1e-13）
 */
void ApplyScalingLimiter(FluidSubCell1d& cell, double epsilon = 1e-13);

/**
 * @brief 对全场所有宏观单元应用缩放限制器
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param block     [in]     网格块信息
 * @param epsilon   [in]     安全下限阈值
 */
void ApplyScalingLimiterAll(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double epsilon = 1e-13
);
