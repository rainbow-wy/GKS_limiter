/**
 * @file interface_1d_limiter.cpp
 * @brief 1D界面通量计算与保正限制实现
 */

#include "interface_1d_limiter.h"
#include <cmath>
#include <algorithm>

//=============================================================================
// 辅助函数：VanLeer限制器
//=============================================================================

static double VanLeerLimiter(double a, double b) {
    if (a * b > 0.0) {
        return 2.0 * a * b / (a + b);
    }
    return 0.0;
}

//=============================================================================
// 宏观界面通量计算
//=============================================================================

void ComputeMacroFluxHigh(
    double* flux,
    Fluid1d* fluids,
    int iface_idx,
    double dt
) {
    // 保存全局变量
    GKS1d_type saved_solver = gks1dsolver;
    
    // 设置为二阶GKS
    gks1dsolver = gks2nd;
    
    // 创建Interface1d结构
    Interface1d iface;
    
    // 调用现有WENO5_AO重构
    // 界面 iface_idx 在 cell[iface_idx-1] 和 cell[iface_idx] 之间
    // 需要从两侧单元分别重构：
    //   - 从左单元 (iface_idx-1) 的右边界获取界面左值
    //   - 从右单元 (iface_idx) 的左边界获取界面右值
    
    // 从左单元向右外推到界面
    Point1d dummy_left;
    WENO5_AO(dummy_left, iface.left, &fluids[iface_idx - 1]);
    // iface.left = 单元 (iface_idx-1) 的右边界值
    
    // 从右单元向左外推到界面
    Point1d dummy_right;
    WENO5_AO(iface.right, dummy_right, &fluids[iface_idx]);
    // iface.right = 单元 (iface_idx) 的左边界值
    
    // 调用g0重构 (Center_collision)
    Center_collision(iface, &fluids[iface_idx - 1]);
    
    // 计算通量
    Flux1d f;
    GKS(f, iface, dt);
    
    // 复制结果
    Copy3(flux, f.f);
    
    // 恢复全局变量
    gks1dsolver = saved_solver;
}

void ComputeMacroFluxLow(
    double* flux,
    Fluid1d* fluids,
    int iface_idx,
    double dt
) {
    // 保存全局变量
    GKS1d_type saved_solver = gks1dsolver;
    
    // 设置为一阶KFVS (保正)
    gks1dsolver = kfvs1st;
    
    // 创建Interface1d结构
    Interface1d iface;
    
    // 一阶重构：直接使用单元平均值
    // 界面 iface_idx 在 cell[iface_idx-1] 和 cell[iface_idx] 之间
    for (int v = 0; v < 3; v++) {
        iface.left.convar[v] = fluids[iface_idx - 1].convar[v];
        iface.right.convar[v] = fluids[iface_idx].convar[v];
        iface.left.der1[v] = 0.0;
        iface.right.der1[v] = 0.0;
    }
    
    // 简单平均作为center值
    for (int v = 0; v < 3; v++) {
        iface.center.convar[v] = 0.5 * (iface.left.convar[v] + iface.right.convar[v]);
        iface.center.der1[v] = 0.0;
    }
    
    // 计算通量
    Flux1d f;
    GKS(f, iface, dt);
    
    // 复制结果
    Copy3(flux, f.f);
    
    // 恢复全局变量
    gks1dsolver = saved_solver;
}

void BlendFlux(
    double* flux_blend,
    const double* flux_high,
    const double* flux_low,
    double alpha_L,
    double alpha_R
) {
    double alpha_bar = 0.5 * (alpha_L + alpha_R);
    double one_minus_alpha = 1.0 - alpha_bar;
    
    for (int v = 0; v < 3; v++) {
        flux_blend[v] = one_minus_alpha * flux_high[v] + alpha_bar * flux_low[v];
    }
}

void TentativeUpdate(
    double* u_L_new,
    double* u_R_new,
    const SubCell1d& sc_L,
    const SubCell1d& sc_R,
    const double* flux,
    double dt
) {
    // 左子单元：通量从右边界流出 (F_left=0, F_right=flux)
    // U = U_old + (1/dx) * (F_left - F_right) = U_old - flux/dx
    double dtdx_L = dt / sc_L.dx;
    for (int v = 0; v < 3; v++) {
        u_L_new[v] = sc_L.convar[v] - dtdx_L * flux[v];  // 右边界通量
    }
    
    // 右子单元：通量从左边界流入 (F_left=flux, F_right=0)
    // U = U_old + (1/dx) * (F_left - F_right) = U_old + flux/dx
    double dtdx_R = dt / sc_R.dx;
    for (int v = 0; v < 3; v++) {
        u_R_new[v] = sc_R.convar[v] + dtdx_R * flux[v];  // 左边界通量
    }
}

//=============================================================================
// 辅助函数：从守恒变量计算压力
//=============================================================================

static double ComputePressure(const double* convar) {
    double rho = convar[0];
    if (rho <= 0.0) return -1.0;  // 无效
    double rho_u = convar[1];
    double E = convar[2];
    double u = rho_u / rho;
    return (Gamma - 1.0) * (E - 0.5 * rho * u * u);
}

//=============================================================================
// 辅助函数：计算单个变量的限制 theta
//=============================================================================

static double ComputeThetaForVariable(
    double val_new,   // 用当前 F_limited 更新后的值
    double val_low,   // 用 F_low 更新后的安全值
    double eps        // 保正阈值
) {
    if (val_new >= eps) {
        // 已满足保正性，不需要限制
        return 1.0;
    }
    
    // val_new < eps，需要限制
    // 公式: theta = |( eps - val_low ) / ( val_low - val_new )|
    double denom = val_low - val_new;
    
    if (std::abs(denom) < 1e-14) {
        // 分母太小，使用安全值
        return 0.0;
    }
    
    double ratio = (eps - val_low) / denom;
    ratio = std::abs(ratio);
    
    return std::min(1.0, std::max(0.0, ratio));
}

//=============================================================================
// 顺序限制策略：先密度，再压力
//=============================================================================

void ComputeMacroFluxFinal(
    double* flux_final,
    const double* flux_high,
    const double* flux_low,
    const SubCell1d& sc_L,
    const SubCell1d& sc_R,
    double alpha_L,
    double alpha_R,
    double dt
) {
    const double eps = POSITIVITY_EPS;
    
    //=========================================================================
    // Step 1: 初始化 F_limited = BlendFlux(F_high, F_low, alpha)
    //=========================================================================
    double F_limited[3];
    BlendFlux(F_limited, flux_high, flux_low, alpha_L, alpha_R);
    
    double lambda_L = 1.0 / sc_L.dx;
    double lambda_R = 1.0 / sc_R.dx;
    
    //=========================================================================
    // Step 2: 顺序限制循环 k = 0 (密度), k = 1 (压力)
    //=========================================================================
    for (int k = 0; k < 2; k++) {
        
        // --- 计算用 F_limited 更新后的状态 ---
        double u_L_new[3], u_R_new[3];
        for (int v = 0; v < 3; v++) {
            u_L_new[v] = sc_L.convar[v] - lambda_L * F_limited[v];  // 左：通量流出
            u_R_new[v] = sc_R.convar[v] + lambda_R * F_limited[v];  // 右：通量流入
        }
        
        // --- 计算用 F_low 更新后的安全状态 ---
        double u_L_low[3], u_R_low[3];
        for (int v = 0; v < 3; v++) {
            u_L_low[v] = sc_L.convar[v] - lambda_L * flux_low[v];
            u_R_low[v] = sc_R.convar[v] + lambda_R * flux_low[v];
        }
        
        // --- 评估物理量 ---
        double val_L_new, val_L_low;
        double val_R_new, val_R_low;
        
        if (k == 0) {
            // 密度检查
            val_L_new = u_L_new[0];
            val_L_low = u_L_low[0];
            val_R_new = u_R_new[0];
            val_R_low = u_R_low[0];
        } else {
            // 压力检查
            val_L_new = ComputePressure(u_L_new);
            val_L_low = ComputePressure(u_L_low);
            val_R_new = ComputePressure(u_R_new);
            val_R_low = ComputePressure(u_R_low);
        }
        
        // --- 计算各子单元的 theta ---
        double theta_L = ComputeThetaForVariable(val_L_new, val_L_low, eps);
        double theta_R = ComputeThetaForVariable(val_R_new, val_R_low, eps);
        
        // --- 取全局最小 theta ---
        double theta_k = std::min(theta_L, theta_R);
        
        // --- 立即更新 F_limited ---
        if (theta_k < 1.0) {
            for (int v = 0; v < 3; v++) {
                F_limited[v] = theta_k * F_limited[v] + (1.0 - theta_k) * flux_low[v];
            }
        }
    }
    
    //=========================================================================
    // Step 3: 输出最终通量
    //=========================================================================
    Copy3(flux_final, F_limited);
}

//=============================================================================
// 子单元内部通量计算
//=============================================================================

void SubcellLinearRecon(
    double* limited_slope,
    const SubCell1d& sc_left,
    const SubCell1d& sc_center,
    const SubCell1d& sc_right
) {
    // 基于 (L, C, R) 三点计算 C 单元的受限斜率
    
    for (int v = 0; v < 3; v++) {
        // 计算两侧原始斜率
        double dx_lc = 0.5 * (sc_left.dx + sc_center.dx);
        double dx_cr = 0.5 * (sc_center.dx + sc_right.dx);
        
        double slope_l = (sc_center.convar[v] - sc_left.convar[v]) / dx_lc;
        double slope_r = (sc_right.convar[v] - sc_center.convar[v]) / dx_cr;
        
        // VanLeer 限制器
        limited_slope[v] = VanLeerLimiter(slope_l, slope_r);
        
        // --- 强激波保护 (Shock Guard) ---
        // 针对 Sod 问题初始的无限梯度，如果左右密度差太大，强制降阶为一阶
        if (v == 0) { // 仅检测密度
            double rho_L = sc_left.convar[0];
            double rho_R = sc_right.convar[0];
            double rho_C = sc_center.convar[0];
            
            double delta_max = std::max(std::abs(rho_L - rho_C), std::abs(rho_R - rho_C));
            double rho_base = rho_C + 1e-10;
            
            // 如果局部跳变超过 20%，强行抹平斜率 (回退到一阶)
            if (delta_max / rho_base > 0.2) {
                limited_slope[0] = 0.0; 
                limited_slope[1] = 0.0; 
                limited_slope[2] = 0.0;
                break; // 动量和能量已清零，跳出循环
            }
        }
    }
}

void ComputeInternalFlux(
    double* flux,
    const SubCell1d& sc_left,
    const SubCell1d& sc_right,
    const double* left_der,
    const double* right_der,
    double dt
) {
    // 保存全局变量
    GKS1d_type saved_solver = gks1dsolver;
    
    // 设置为二阶KFVS
    gks1dsolver = kfvs2nd;
    
    // 创建Interface1d结构
    Interface1d iface;
    
    // 填充左右值和斜率
    for (int v = 0; v < 3; v++) {
        // 左值：sc_left外推到右边界
        iface.left.convar[v] = sc_left.convar[v] + 0.5 * sc_left.dx * left_der[v];
        // 右值：sc_right外推到左边界
        iface.right.convar[v] = sc_right.convar[v] - 0.5 * sc_right.dx * right_der[v];
        
        // 斜率
        iface.left.der1[v] = left_der[v];
        iface.right.der1[v] = right_der[v];
    }
    
    // center值 (简单平均)
    for (int v = 0; v < 3; v++) {
        iface.center.convar[v] = 0.5 * (iface.left.convar[v] + iface.right.convar[v]);
        iface.center.der1[v] = 0.5 * (iface.left.der1[v] + iface.right.der1[v]);
    }
    
    // 计算通量
    Flux1d f;
    GKS(f, iface, dt);
    
    // 复制结果
    Copy3(flux, f.f);
    
    // 恢复全局变量
    gks1dsolver = saved_solver;
}

//=============================================================================
// 全场通量计算
//=============================================================================

void ComputeAllMacroFluxes(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block,
    double dt
) {
    // 遍历所有宏观界面
    // 界面 i 在 cell[i-1] 和 cell[i] 之间
    for (int i = block.ghost; i <= block.nx - block.ghost; i++) {
        // 计算高阶通量
        double flux_high[3];
        ComputeMacroFluxHigh(flux_high, fluids, i, dt);
        
        // 计算低阶通量
        double flux_low[3];
        ComputeMacroFluxLow(flux_low, fluids, i, dt);
        
        // 获取界面两侧的子单元
        FluidSubCell1d& cell_L = fluids_sc[i - 1];
        FluidSubCell1d& cell_R = fluids_sc[i];
        
        // 计算最终通量
        double flux_final[3];
        ComputeMacroFluxFinal(
            flux_final,
            flux_high, flux_low,
            cell_L.sub[N_SUBCELL - 1],  // 左单元最右子单元
            cell_R.sub[0],              // 右单元最左子单元
            cell_L.alpha, cell_R.alpha,
            dt
        );
        
        // 存储通量
        Copy3(cell_L.flux_right, flux_final);
        Copy3(cell_R.flux_left, flux_final);
        
        // 存储高阶/低阶通量 (调试用)
        Copy3(cell_L.flux_high_right, flux_high);
        Copy3(cell_L.flux_low_right, flux_low);
    }
}

void ComputeAllInternalFluxes(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double dt
) {
    // 遍历所有宏观单元
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        FluidSubCell1d& cell = fluids_sc[i];
        
        // 计算每个内部子界面的通量 (共 N_SUBCELL-1 = 4 个)
        for (int j = 0; j < N_SUBCELL - 1; j++) {
            // 界面 j 在 sub[j] 和 sub[j+1] 之间
            
            // 重新计算受限斜率：
            // 界面 j 在 sub[j] 和 sub[j+1] 之间
            double slope_L[3], slope_R[3];
            
            // 左侧 sub[j] 的斜率：使用 (j-1, j, j+1)
            const SubCell1d& sc_L_prev = (j > 0) ? cell.sub[j - 1] : cell.sub[j];
            const SubCell1d& sc_L_curr = cell.sub[j];
            const SubCell1d& sc_L_next = cell.sub[j + 1];
            SubcellLinearRecon(slope_L, sc_L_prev, sc_L_curr, sc_L_next);
            
            // 右侧 sub[j+1] 的斜率：使用 (j, j+1, j+2)
            const SubCell1d& sc_R_prev = cell.sub[j];
            const SubCell1d& sc_R_curr = cell.sub[j + 1];
            const SubCell1d& sc_R_next = (j + 2 < N_SUBCELL) ? cell.sub[j + 2] : cell.sub[j + 1];
            SubcellLinearRecon(slope_R, sc_R_prev, sc_R_curr, sc_R_next);
            
            // 计算通量 (传入限制后的斜率)
            ComputeInternalFlux(
                cell.flux_internal[j],
                cell.sub[j], cell.sub[j + 1],
                slope_L, slope_R,
                dt
            );
        }
    }
}

//=============================================================================
// 缩放限制器 (Scaling Limiter) - Zhang & Shu, 2010 算法实现
//=============================================================================

/**
 * @brief 从守恒变量计算压力（内部辅助函数）
 * 
 * @param convar 守恒变量 [ρ, m, E]
 * @return 压力 p = (γ-1)(E - 0.5*m²/ρ)
 */
static double ComputePressureFromConvar(const double* convar) {
    double rho = convar[0];
    if (rho <= 0.0) return -1.0;  // 无效密度
    double m = convar[1];
    double E = convar[2];
    return (Gamma - 1.0) * (E - 0.5 * m * m / rho);
}

void ApplyScalingLimiter(FluidSubCell1d& cell, double epsilon) {
    
    //=========================================================================
    // 步骤 0：计算宏观平均值 U_mean（作为安全锚点）
    // 使用 GL 加权平均：U_mean = Σ (weight[j]/2) * U_sub[j]
    //=========================================================================
    double U_mean[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < N_SUBCELL; j++) {
        double w = cell.sub[j].weight / 2.0;
        for (int v = 0; v < 3; v++) {
            U_mean[v] += w * cell.sub[j].convar[v];
        }
    }
    
    double rho_avg = U_mean[0];
    double m_avg   = U_mean[1];
    double E_avg   = U_mean[2];
    
    //=========================================================================
    // 步骤 1：密度限制（线性修正）
    //=========================================================================
    
    // 1.1 找到最小密度
    double rho_min = cell.sub[0].convar[0];
    for (int j = 1; j < N_SUBCELL; j++) {
        rho_min = std::min(rho_min, cell.sub[j].convar[0]);
    }
    
    // 1.2 计算密度缩放系数 θ_ρ
    double theta_rho = 1.0;
    if (rho_min < epsilon) {
        // θ_ρ = (ρ_avg - ε) / (ρ_avg - ρ_min)
        double denom = rho_avg - rho_min;
        if (std::abs(denom) > 1e-15) {
            theta_rho = (rho_avg - epsilon) / denom;
            // 确保 θ_ρ ∈ [0, 1]
            theta_rho = std::max(0.0, std::min(1.0, theta_rho));
        } else {
            // 分母太小，回退到平均值
            theta_rho = 0.0;
        }
    }
    
    // 1.3 执行密度缩放：U_sub[j] = U_mean + θ_ρ * (U_sub[j] - U_mean)
    if (theta_rho < 1.0) {
        for (int j = 0; j < N_SUBCELL; j++) {
            for (int v = 0; v < 3; v++) {
                cell.sub[j].convar[v] = U_mean[v] + theta_rho * (cell.sub[j].convar[v] - U_mean[v]);
            }
        }
    }
    
    //=========================================================================
    // 步骤 2：压力限制（二次方程修正）
    // 目标：确保 P = (γ-1)(E - 0.5*m²/ρ) >= ε
    // 等价于：2ρE - m² - 2ρε/(γ-1) >= 0
    //=========================================================================
    
    double theta_p = 1.0;  // 全局压力缩放系数
    
    double gamma_m1 = Gamma - 1.0;  // γ - 1
    
    for (int j = 0; j < N_SUBCELL; j++) {
        // 获取当前子单元状态
        double rho_j = cell.sub[j].convar[0];
        double m_j   = cell.sub[j].convar[1];
        double E_j   = cell.sub[j].convar[2];
        
        // 计算当前压力
        double P_j = (rho_j > epsilon) ? 
                     (gamma_m1 * (E_j - 0.5 * m_j * m_j / rho_j)) : -1.0;
        
        // 如果压力已满足条件，跳过
        if (P_j >= epsilon) {
            continue;
        }
        
        // 需要求解：找到 θ 使得 P(U_mean + θ*(U_j - U_mean)) = ε
        // 偏差向量 D = U_j - U_mean
        double rho_d = rho_j - rho_avg;
        double m_d   = m_j - m_avg;
        double E_d   = E_j - E_avg;
        
        // 二次方程系数：Aθ² + Bθ + C = 0
        // A = 2*ρ_d*E_d - m_d²
        // B = 2*(ρ_avg*E_d + ρ_d*E_avg) - 2*m_avg*m_d
        // C = 2*ρ_avg*E_avg - m_avg² - 2*ρ_avg*ε/(γ-1)
        double A = 2.0 * rho_d * E_d - m_d * m_d;
        double B = 2.0 * (rho_avg * E_d + rho_d * E_avg) - 2.0 * m_avg * m_d;
        double C = 2.0 * rho_avg * E_avg - m_avg * m_avg - 2.0 * rho_avg * epsilon / gamma_m1;
        
        // 求解二次方程
        double theta_local = 0.0;  // 默认回退到平均值
        
        if (std::abs(A) < 1e-15) {
            // A ≈ 0，退化为线性方程：Bθ + C = 0
            if (std::abs(B) > 1e-15) {
                theta_local = -C / B;
            }
        } else {
            // 标准二次方程求解
            double discriminant = B * B - 4.0 * A * C;
            
            if (discriminant >= 0.0) {
                double sqrt_disc = std::sqrt(discriminant);
                // 两个根：(-B ± sqrt(Δ)) / (2A)
                double root1 = (-B - sqrt_disc) / (2.0 * A);
                double root2 = (-B + sqrt_disc) / (2.0 * A);
                
                // 选择在 [0, 1] 区间内、且靠近 1 的根
                // 注意：θ=1 对应当前状态（P<ε），θ=0 对应平均值（P_avg>ε）
                // 我们需要找到使 P=ε 的正根
                
                // 判断哪个根是有效的（在 [0, 1] 区间内）
                bool valid1 = (root1 >= 0.0 && root1 <= 1.0);
                bool valid2 = (root2 >= 0.0 && root2 <= 1.0);
                
                if (valid1 && valid2) {
                    // 两个根都有效，选择较小的（更安全）
                    theta_local = std::min(root1, root2);
                } else if (valid1) {
                    theta_local = root1;
                } else if (valid2) {
                    theta_local = root2;
                } else {
                    // 没有有效根在 [0,1] 内
                    // 这不应该发生（因为平均值是安全的），回退到 0
                    theta_local = 0.0;
                }
            } else {
                // 判别式为负，无实根，回退到平均值
                theta_local = 0.0;
            }
        }
        
        // 确保 θ_local ∈ [0, 1]
        theta_local = std::max(0.0, std::min(1.0, theta_local));
        
        // 更新全局最小系数
        theta_p = std::min(theta_p, theta_local);
    }
    
    // 2.2 执行压力缩放：U_sub[j] = U_mean + θ_p * (U_sub[j] - U_mean)
    if (theta_p < 1.0) {
        for (int j = 0; j < N_SUBCELL; j++) {
            for (int v = 0; v < 3; v++) {
                cell.sub[j].convar[v] = U_mean[v] + theta_p * (cell.sub[j].convar[v] - U_mean[v]);
            }
        }
    }
}

void ApplyScalingLimiterAll(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double epsilon
) {
    // 遍历所有内部宏观单元（不含 ghost 单元）
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        ApplyScalingLimiter(fluids_sc[i], epsilon);
    }
}

