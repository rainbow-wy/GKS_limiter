/**
 * @file subcell_1d_tools.h
 * @brief 1D子单元工具函数接口
 * 
 * 包含Legendre投影和光滑度因子计算
 */

#pragma once
#include "subcell_1d_data.h"

//=============================================================================
// 光滑度因子计算参数 (硬编码)
//=============================================================================

/// 阈值参数 a
const double ALPHA_PARAM_A = 0.5;

/// 阈值参数 c  
const double ALPHA_PARAM_C = 1.8;

/// Logistic函数锐度系数 S
const double ALPHA_PARAM_S = 9.21024;

/// 截断阈值
const double ALPHA_MIN = 0.001;

//=============================================================================
// Legendre多项式
//=============================================================================

/**
 * @brief 计算Legendre多项式 P_n(x)
 * 
 * 使用递推公式: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
 * 
 * @param n 多项式阶数 (0, 1, 2, ...)
 * @param x 自变量 ∈ [-1, 1]
 * @return P_n(x) 的值
 */
double Legendre(int n, double x);

//=============================================================================
// 光滑度因子计算
//=============================================================================

/**
 * @brief 计算单个宏观单元的Legendre系数
 * 
 * 公式: q̂_j = Σ_{k=0}^{4} q(ξ_k) * L_j(2*ξ_k - 1) * w_k
 * 其中 q = 密度
 * 
 * @param cell [in/out] 子单元扩展结构（legendre_coeff被更新）
 */
void ComputeLegendreCoeffs(FluidSubCell1d& cell);

/**
 * @brief 计算模态能量比 E
 * 
 * 公式: E = max( q̂_3² / Σ(j=0..3)q̂_j², q̂_4² / Σ(j=0..4)q̂_j² )
 * 
 * @param cell [in] 子单元扩展结构
 * @return 能量比 E
 */
double ComputeModalEnergy(const FluidSubCell1d& cell);

/**
 * @brief 计算阈值 T(N)
 * 
 * 公式: T = a * 10^(-c * (N+1)^0.25)
 * 其中 a=0.5, c=1.8, N+1=5
 * 
 * @return 阈值 T
 */
double ComputeThreshold();

/**
 * @brief Logistic映射计算原始α
 * 
 * 公式: α_raw = 1 / (1 + exp(-S/T * (E - T)))
 * 
 * @param E 能量比
 * @param T 阈值
 * @return 原始α值
 */
double LogisticMapping(double E, double T);

/**
 * @brief 截断α值
 * 
 * 规则:
 * - if α < 0.001: α = 0
 * - elif α > 0.999: α = 1
 * - else: α不变
 * 
 * @param alpha [in] 原始α值
 * @return 截断后的α值
 */
double ClipAlpha(double alpha);

/**
 * @brief 计算单个宏观单元的光滑度因子
 * 
 * 完整流程:
 * 1. Legendre投影
 * 2. 计算能量比E
 * 3. Logistic映射
 * 4. 截断
 * 
 * @param cell [in/out] 子单元扩展结构（alpha被更新）
 */
void ComputeSmoothnessIndicator(FluidSubCell1d& cell);

/**
 * @brief 平滑光滑度因子 (考虑相邻单元)
 * 
 * 公式: α_final = max(α, 0.5*α_left, 0.5*α_right)
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param block     [in]     网格块信息
 */
void SmoothAlpha(FluidSubCell1d* fluids_sc, Block1d& block);

/**
 * @brief 计算全场光滑度因子
 * 
 * 包含计算+平滑两个步骤
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param block     [in]     网格块信息
 */
void ComputeAllSmoothnessIndicators(FluidSubCell1d* fluids_sc, Block1d& block);
