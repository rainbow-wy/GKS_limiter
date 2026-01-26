/**
 * @file gks_subcell_1d_solver.h
 * @brief 1D子单元GKS求解器主接口
 * 
 * 完整求解流程：初始化 -> 重构 -> 光滑度 -> 通量 -> 演化
 */

#pragma once
#include "subcell_1d_data.h"
#include "subcell_1d_tools.h"
#include "interface_1d_limiter.h"

//=============================================================================
// 求解器初始化与销毁
//=============================================================================

/**
 * @brief 初始化1D子单元求解器
 * 
 * @param block  [in]  网格块信息
 * @param fluids [in]  宏观单元数组
 * @return 分配并初始化的FluidSubCell1d数组
 */
FluidSubCell1d* InitSubCell1dSolver(
    Block1d& block,
    Fluid1d* fluids
);

/**
 * @brief 销毁1D子单元求解器
 * 
 * @param fluids_sc [in] 子单元数组
 */
void FinalizeSubCell1dSolver(FluidSubCell1d* fluids_sc);

//=============================================================================
// 演化函数
//=============================================================================

/**
 * @brief 单时间步演化
 * 
 * 完整流程:
 * 1. 保存old状态
 * 2. 计算光滑度因子α
 * 3. 计算宏观界面通量
 * 4. 计算子单元内部通量
 * 5. 更新子单元 (非均匀体积)
 * 6. 投影回宏观单元
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param fluids    [in/out] 宏观单元数组
 * @param block     [in]     网格块信息
 */
void SubCell1dSolverStep(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
);

/**
 * @brief 更新所有子单元
 * 
 * U^{n+1} = U^n - dt/dx_sub * (F_right - F_left)
 * 注意: dx_sub 是基于GL权重的非均匀网格
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param block     [in]     网格块信息
 * @param dt        [in]     时间步长
 */
void UpdateAllSubCells1d(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double dt
);

/**
 * @brief 全场投影到宏观单元
 * 
 * @param fluids_sc [in]     子单元数组
 * @param fluids    [out]    宏观单元数组
 * @param block     [in]     网格块信息
 */
void ProjectAllToMacro1d(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
);

//=============================================================================
// 边界条件
//=============================================================================

/**
 * @brief 将宏观单元边界条件传递到子单元
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param fluids    [in]     宏观单元数组 (已设置好边界)
 * @param block     [in]     网格块信息
 */
void ApplySubCell1dBC(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
);

//=============================================================================
// 输出函数
//=============================================================================

/**
 * @brief 输出子单元数据到.plt文件
 * 
 * 格式: variables = x, density, u, pressure, temperature, entropy, Ma, alpha
 * 
 * @param fluids_sc [in] 子单元数组
 * @param block     [in] 网格块信息
 */
void OutputSubCell1d(FluidSubCell1d* fluids_sc, Block1d& block);

/**
 * @brief 输出通量调试信息到CSV文件
 * 
 * 格式: macro_idx, sub_idx, x, rho, rho_old, flux_left, flux_right, alpha
 * 
 * @param fluids_sc [in] 子单元数组
 * @param block     [in] 网格块信息
 */
void OutputFluxDebug(FluidSubCell1d* fluids_sc, Block1d& block);

/**
 * @brief 打印单个宏观单元的子单元信息 (调试用)
 * 
 * @param cell    [in] 子单元扩展结构
 * @param verbose [in] 是否输出详细信息
 */
void PrintSubCell1dInfo(const FluidSubCell1d& cell, bool verbose = false);
