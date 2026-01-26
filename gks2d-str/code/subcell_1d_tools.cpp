/**
 * @file subcell_1d_tools.cpp
 * @brief 1D子单元工具函数实现
 * 
 * Legendre投影和光滑度因子计算
 */

#include "subcell_1d_tools.h"
#include <cmath>
#include <algorithm>

//=============================================================================
// Legendre多项式
//=============================================================================

double Legendre(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    double Pnm1 = 1.0;  // P_{n-1}
    double Pn = x;      // P_n
    double Pnp1;        // P_{n+1}
    
    for (int k = 1; k < n; k++) {
        // 递推公式: (k+1)P_{k+1}(x) = (2k+1)xP_k(x) - kP_{k-1}(x)
        Pnp1 = ((2.0 * k + 1.0) * x * Pn - k * Pnm1) / (k + 1.0);
        Pnm1 = Pn;
        Pn = Pnp1;
    }
    
    return Pn;
}

//=============================================================================
// 光滑度因子计算
//=============================================================================

void ComputeLegendreCoeffs(FluidSubCell1d& cell) {
    // 公式: q̂_j = Σ_{k=0}^{4} q(ξ_k) * L_j(2*ξ_k - 1) * w_k
    // 其中 q = 子单元密度, ξ_k 为GL点位置
    
    // 清零
    for (int j = 0; j <= N_SUBCELL; j++) {
        cell.legendre_coeff[j] = 0.0;
    }
    
    // 遍历每个子单元（即每个GL积分点）
    for (int k = 0; k < N_SUBCELL; k++) {
        double xi_k = GL_POINTS_5[k];       // 参考坐标 ξ ∈ [-1, 1]
        double w_k = GL_WEIGHTS_5[k];       // GL权重
        double q_k = cell.sub[k].convar[0]; // 密度作为指示变量
        
        // 坐标映射: ξ ∈ [-1, 1] 直接使用（无需额外变换）
        // Legendre多项式定义在 [-1, 1] 上
        
        // 对每阶Legendre系数累加
        for (int j = 0; j <= N_SUBCELL - 1; j++) {  // j = 0, 1, 2, 3, 4
            double L_j = Legendre(j, xi_k);
            cell.legendre_coeff[j] += q_k * L_j * w_k;
        }
    }
    
    // j = 0, 1, 2, 3, 4 are covered in the loop above
}

double ComputeModalEnergy(const FluidSubCell1d& cell) {
    // 公式: E = max( q̂_3² / Σ(j=0..3)q̂_j², q̂_4² / Σ(j=0..4)q̂_j² )
    
    double eps = 1e-14;  // 防止除零
    
    // 计算 q̂_3² / Σ(j=0..3)q̂_j²
    double sum_03 = eps;
    for (int j = 0; j <= 3; j++) {
        sum_03 += cell.legendre_coeff[j] * cell.legendre_coeff[j];
    }
    double E1 = cell.legendre_coeff[3] * cell.legendre_coeff[3] / sum_03;
    
    // 计算 q̂_4² / Σ(j=0..4)q̂_j²
    double sum_04 = eps;
    for (int j = 0; j <= 4; j++) {
        sum_04 += cell.legendre_coeff[j] * cell.legendre_coeff[j];
    }
    double E2 = cell.legendre_coeff[4] * cell.legendre_coeff[4] / sum_04;
    
    return std::max(E1, E2);
}

double ComputeThreshold() {
    // 公式: T = a * 10^(-c * (N+1)^0.25)
    // 其中 a = 0.5, c = 1.8, N+1 = 5
    
    double exponent = -ALPHA_PARAM_C * std::pow(5.0, 0.25);
    return ALPHA_PARAM_A * std::pow(10.0, exponent);
}

double LogisticMapping(double E, double T) {
    // 公式: α_raw = 1 / (1 + exp(-S/T * (E - T)))
    
    double arg = -ALPHA_PARAM_S / T * (E - T);
    
    // 防止溢出
    if (arg > 20.0) return 0.0;
    if (arg < -20.0) return 1.0;
    
    return 1.0 / (1.0 + std::exp(arg));
}

double ClipAlpha(double alpha) {
    // 截断规则:
    // - if α < 0.001: α = 0
    // - elif α > 0.999: α = 1
    // - else: α不变
    
    if (alpha < ALPHA_MIN) {
        return 0.0;
    } else if (alpha > 1.0 - ALPHA_MIN) {
        return 1.0;
    } else {
        return alpha;
    }
}

void ComputeSmoothnessIndicator(FluidSubCell1d& cell) {
    // Step 1: Legendre投影
    ComputeLegendreCoeffs(cell);
    
    // Step 2: 计算能量比E
    double E = ComputeModalEnergy(cell);
    
    // Step 3: 计算阈值T
    double T = ComputeThreshold();
    
    // Step 4: Logistic映射
    cell.alpha_raw = LogisticMapping(E, T);
    
    // Step 5: 截断
    cell.alpha = ClipAlpha(cell.alpha_raw);
}

void SmoothAlpha(FluidSubCell1d* fluids_sc, Block1d& block) {
    // 公式: α_final = max(α, 0.5*α_left, 0.5*α_right)
    
    // 需要临时存储原始α值
    double* alpha_orig = new double[block.nx];
    for (int i = 0; i < block.nx; i++) {
        alpha_orig[i] = fluids_sc[i].alpha;
    }
    
    // 对内部单元进行平滑
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        double alpha_self = alpha_orig[i];
        double alpha_left = (i > 0) ? alpha_orig[i - 1] : alpha_orig[i];
        double alpha_right = (i < block.nx - 1) ? alpha_orig[i + 1] : alpha_orig[i];
        
        fluids_sc[i].alpha = std::max({
            alpha_self,
            0.5 * alpha_left,
            0.5 * alpha_right
        });
    }
    
    delete[] alpha_orig;
}

void ComputeAllSmoothnessIndicators(FluidSubCell1d* fluids_sc, Block1d& block) {
    // Step 1: 计算每个单元的α
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        ComputeSmoothnessIndicator(fluids_sc[i]);
    }
    
    // Step 2: 平滑
    SmoothAlpha(fluids_sc, block);
}
