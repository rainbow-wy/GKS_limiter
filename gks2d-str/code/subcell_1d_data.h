/**
 * @file subcell_1d_data.h
 * @brief 1D子单元数据结构定义
 * 
 * 基于5点Gauss-Legendre积分点划分子单元
 * 
 * @note 遵守Zero-Modification原则，不修改任何现有源文件
 */

#pragma once
#include "fvm_gks1d.h"

//=============================================================================
// 常量定义
//=============================================================================

/// 子单元数量 (N+1 = 5, 对应5点GL积分)
const int N_SUBCELL = 5;

/// 5点 Gauss-Legendre 积分权重 (硬编码)
const double GL_WEIGHTS_5[5] = {
    0.236926885056189,   // ξ = -0.906179845938664
    0.478628670499366,   // ξ = -0.538469310105683
    0.568888888888889,   // ξ =  0.0
    0.478628670499366,   // ξ =  0.538469310105683
    0.236926885056189    // ξ =  0.906179845938664
};

/// 5点 Gauss-Legendre 积分点位置 (参考坐标 ξ ∈ [-1, 1])
const double GL_POINTS_5[5] = {
    -0.906179845938664,
    -0.538469310105683,
     0.0,
     0.538469310105683,
     0.906179845938664
};

/// 保正容差
const double POSITIVITY_EPS = 1e-10;

//=============================================================================
// 1D子单元数据结构
//=============================================================================

/**
 * @brief 单个子单元数据
 */
struct SubCell1d {
    double convar[3];       ///< 守恒变量 [ρ, ρu, E]
    double convar_old[3];   ///< 旧时刻守恒变量
    double primvar[3];      ///< 原始变量 [ρ, u, p]
    
    double dx;              ///< 子单元宽度 = macro_dx * weight / 2
    double center_x;        ///< 子单元中心坐标
    double weight;          ///< GL权重
    
    int local_idx;          ///< 在宏观单元内的索引 (0 ~ N_SUBCELL-1)
};

/**
 * @brief 扩展的宏观单元（包含子单元）
 */
struct FluidSubCell1d {
    Fluid1d* parent;        ///< 指向父宏观单元
    int macro_idx;          ///< 宏观单元在Block中的索引
    
    SubCell1d sub[N_SUBCELL];   ///< 子单元数组
    
    //-------------------------------------------------------------------------
    // 光滑度因子
    //-------------------------------------------------------------------------
    double legendre_coeff[N_SUBCELL + 1];  ///< Legendre系数 q̂_0 ~ q̂_4
    double alpha;           ///< 光滑度因子 α_e ∈ [0, 1]
    double alpha_raw;       ///< 原始α (截断前)
    
    //-------------------------------------------------------------------------
    // 通量存储
    //-------------------------------------------------------------------------
    /// 左宏观边界最终通量
    double flux_left[3];
    /// 右宏观边界最终通量  
    double flux_right[3];
    /// 左宏观边界高阶通量
    double flux_high_left[3];
    /// 左宏观边界低阶通量
    double flux_low_left[3];
    /// 右宏观边界高阶通量
    double flux_high_right[3];
    /// 右宏观边界低阶通量
    double flux_low_right[3];
    /// 内部子界面通量 (共 N_SUBCELL-1 = 4 个)
    double flux_internal[N_SUBCELL - 1][3];
};

//=============================================================================
// 初始化函数
//=============================================================================

/**
 * @brief 分配子单元数组
 * 
 * @param block 1D网格块信息
 * @return 分配的FluidSubCell1d数组
 */
FluidSubCell1d* AllocateSubCell1dArray(Block1d& block);

/**
 * @brief 释放子单元数组
 * 
 * @param fluids_sc 子单元数组指针
 */
void FreeSubCell1dArray(FluidSubCell1d* fluids_sc);

/**
 * @brief 初始化单个宏观单元的子单元几何
 * 
 * 核心公式: sub_dx = macro_dx * weight[i] / 2.0
 * (除以2是因为参考区间 ξ∈[-1,1] 长度为2)
 * 
 * @param cell      [out] 子单元扩展结构
 * @param parent    [in]  父宏观单元
 * @param macro_idx [in]  宏观单元索引
 */
void InitSubCell1dGeometry(
    FluidSubCell1d& cell,
    Fluid1d* parent,
    int macro_idx
);

/**
 * @brief 初始化全场子单元
 * 
 * @param fluids_sc [out] 子单元扩展数组
 * @param fluids    [in]  宏观单元数组
 * @param block     [in]  网格块信息
 */
void InitAllSubCells1d(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
);

//=============================================================================
// 数据操作函数
//=============================================================================

/**
 * @brief 从宏观单元分发初始值到子单元
 * 
 * @param cell [in/out] 子单元扩展结构
 */
void DistributeToSubCells1d(FluidSubCell1d& cell);

/**
 * @brief 将子单元投影回宏观单元 (GL加权平均)
 * 
 * @param cell [in/out] 子单元扩展结构
 */
void ProjectToMacro1d(FluidSubCell1d& cell);

/**
 * @brief 保存子单元状态到old
 * 
 * @param cell [in/out] 子单元扩展结构
 */
void CopySubCell1dNewToOld(FluidSubCell1d& cell);

/**
 * @brief 全场保存状态
 * 
 * @param fluids_sc [in/out] 子单元数组
 * @param block     [in]     网格块信息
 */
void CopyAllSubCells1dNewToOld(FluidSubCell1d* fluids_sc, Block1d& block);

/**
 * @brief 计算子单元原始变量
 * 
 * @param cell [in/out] 子单元扩展结构
 */
void SubCell1dConvarToPrimvar(FluidSubCell1d& cell);

/**
 * @brief 检查子单元保正性
 * 
 * @param convar [in] 守恒变量
 * @param eps    [in] 容差
 * @return true 如果 ρ > eps 且 p > eps
 */
bool IsAdmissible1d(const double* convar, double eps = POSITIVITY_EPS);

//=============================================================================
// 辅助工具
//=============================================================================

/// 数组复制
inline void Copy3(double* dest, const double* src) {
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
}

/// 数组清零
inline void Zero3(double* arr) {
    arr[0] = 0.0;
    arr[1] = 0.0;
    arr[2] = 0.0;
}
