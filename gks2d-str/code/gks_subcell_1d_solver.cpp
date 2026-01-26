/**
 * @file gks_subcell_1d_solver.cpp
 * @brief 1D子单元GKS求解器实现
 */

#include "gks_subcell_1d_solver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

//=============================================================================
// 求解器初始化与销毁
//=============================================================================

FluidSubCell1d* InitSubCell1dSolver(
    Block1d& block,
    Fluid1d* fluids
) {
    cout << "============================================" << endl;
    cout << "Initializing 1D Sub-cell GKS Solver" << endl;
    cout << "  Macro cells: " << block.nodex << endl;
    cout << "  Ghost cells: " << block.ghost << endl;
    cout << "  Sub-cells per macro: " << N_SUBCELL << endl;
    cout << "============================================" << endl;
    
    // 分配数组
    FluidSubCell1d* fluids_sc = AllocateSubCell1dArray(block);
    
    // 初始化几何
    InitAllSubCells1d(fluids_sc, fluids, block);
    
    // 打印阈值信息
    double T = ComputeThreshold();
    cout << "  Threshold T = " << T << endl;
    cout << "  S = " << ALPHA_PARAM_S << endl;
    
    cout << "Sub-cell Solver initialized successfully." << endl;
    
    return fluids_sc;
}

void FinalizeSubCell1dSolver(FluidSubCell1d* fluids_sc) {
    FreeSubCell1dArray(fluids_sc);
    cout << "1D Sub-cell Solver finalized." << endl;
}

//=============================================================================
// 演化函数
//=============================================================================

void SubCell1dSolverStep(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
) {
    double dt = block.dt;
    
    //=========================================================================
    // Step 0: 保存old状态
    //=========================================================================
    CopyAllSubCells1dNewToOld(fluids_sc, block);
    
    //=========================================================================
    // Step 1: 计算光滑度因子α
    //=========================================================================
    ComputeAllSmoothnessIndicators(fluids_sc, block);
    
    //=========================================================================
    // Step 2: 计算宏观界面通量 (高阶GKS + 低阶KFVS + 混合 + 保正修正)
    //=========================================================================
    ComputeAllMacroFluxes(fluids_sc, fluids, block, dt);
    
    //=========================================================================
    // Step 3: 计算子单元内部通量 (二阶KFVS)
    //=========================================================================
    ComputeAllInternalFluxes(fluids_sc, block, dt);
    
    //=========================================================================
    // Step 4: 更新子单元 (非均匀体积)
    //=========================================================================
    UpdateAllSubCells1d(fluids_sc, block, dt);
    
    //=========================================================================
    // Step 5: 缩放限制器 (Zhang & Shu 2010) - 算法最后一道保正保险
    // 确保每个子单元的密度和压力都为正，同时保持宏观平均值不变
    //=========================================================================
    ApplyScalingLimiterAll(fluids_sc, block);
    
    //=========================================================================
    // Step 6: 投影回宏观单元
    //=========================================================================
    ProjectAllToMacro1d(fluids_sc, fluids, block);
}

void UpdateAllSubCells1d(
    FluidSubCell1d* fluids_sc,
    Block1d& block,
    double dt
) {
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        FluidSubCell1d& cell = fluids_sc[i];
        
        for (int j = 0; j < N_SUBCELL; j++) {
            SubCell1d& sc = cell.sub[j];
            
            // 确定左右通量
            double* F_left;
            double* F_right;
            
            if (j == 0) {
                // 最左子单元：左边界是宏观边界
                F_left = cell.flux_left;
            } else {
                // 内部：使用内部子界面通量
                F_left = cell.flux_internal[j - 1];
            }
            
            if (j == N_SUBCELL - 1) {
                // 最右子单元：右边界是宏观边界
                F_right = cell.flux_right;
            } else {
                // 内部：使用内部子界面通量
                F_right = cell.flux_internal[j];
            }
            
            double dtdx = 1.0 / sc.dx;
            for (int v = 0; v < 3; v++) {
                sc.convar[v] = sc.convar_old[v] - dtdx * (F_right[v] - F_left[v]);
            }
        }
    }
}

void ProjectAllToMacro1d(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
) {
    for (int i = block.ghost; i < block.nx - block.ghost; i++) {
        ProjectToMacro1d(fluids_sc[i]);
    }
}

//=============================================================================
// 边界条件
//=============================================================================

void ApplySubCell1dBC(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
) {
    // 将宏观单元边界值传递到ghost区域的子单元
    
    // 左边界 ghost cells
    for (int i = 0; i < block.ghost; i++) {
        DistributeToSubCells1d(fluids_sc[i]);
    }
    
    // 右边界 ghost cells
    for (int i = block.nx - block.ghost; i < block.nx; i++) {
        DistributeToSubCells1d(fluids_sc[i]);
    }
}

//=============================================================================
// 输出函数
//=============================================================================

void OutputSubCell1d(FluidSubCell1d* fluids_sc, Block1d& block) {
    // 生成文件名 (使用绝对路径，与output.cpp一致)
    stringstream name;
    name << "C:/code/gks2d-str/code/result/SubCell1D-" << setfill('0') << setw(5) << block.step << ".plt";
    string s;
    name >> s;
    
    cout << "Output sub-cell file: " << s << endl;
    
    ofstream out(s);
    
    // 文件头
    out << "variables = x, density, u, pressure, temperature, entropy, Ma, alpha" << endl;
    out << "zone i = " << block.nodex * N_SUBCELL << ", F=POINT" << endl;
    
    // 计算原始变量
    for (int i = block.ghost; i < block.ghost + block.nodex; i++) {
        SubCell1dConvarToPrimvar(fluids_sc[i]);
    }
    
    // 输出每个子单元
    for (int i = block.ghost; i < block.ghost + block.nodex; i++) {
        FluidSubCell1d& cell = fluids_sc[i];
        
        for (int j = 0; j < N_SUBCELL; j++) {
            SubCell1d& sc = cell.sub[j];
            
            double rho = sc.primvar[0];
            double u = sc.primvar[1];
            double p = sc.primvar[2];
            
            double T_temp = Temperature(rho, p);
            double s_entropy = entropy(rho, p);
            double c = Soundspeed(rho, p);
            double Ma = (c > 0) ? abs(u) / c : 0.0;
            
            out << setprecision(10)
                << sc.center_x << " "
                << rho << " "
                << u << " "
                << p << " "
                << T_temp << " "
                << s_entropy << " "
                << Ma << " "
                << cell.alpha << endl;
        }
    }
    
    out.close();
}

void OutputFluxDebug(FluidSubCell1d* fluids_sc, Block1d& block) {
    // 生成文件名
    stringstream name;
    name << "C:/code/gks2d-str/code/result/FluxDebug-" << setfill('0') << setw(5) << block.step << ".csv";
    string s;
    name >> s;
    
    cout << "Output flux debug file: " << s << endl;
    
    ofstream out(s);
    
    // CSV 文件头 (添加flux_high_right用于调试)
    out << "macro_idx,sub_idx,x,rho,rho_old,"
        << "flux_left_0,flux_left_1,flux_left_2,"
        << "flux_right_0,flux_right_1,flux_right_2,"
        << "flux_high_right_0,flux_high_right_1,flux_high_right_2,"
        << "flux_low_right_0,flux_low_right_1,flux_low_right_2,"
        << "alpha" << endl;
    
    // 输出每个子单元及其通量
    for (int i = block.ghost; i < block.ghost + block.nodex; i++) {
        FluidSubCell1d& cell = fluids_sc[i];
        
        for (int j = 0; j < N_SUBCELL; j++) {
            SubCell1d& sc = cell.sub[j];
            
            // 确定该子单元的左右通量
            double* F_left;
            double* F_right;
            
            if (j == 0) {
                F_left = cell.flux_left;
            } else {
                F_left = cell.flux_internal[j - 1];
            }
            
            if (j == N_SUBCELL - 1) {
                F_right = cell.flux_right;
            } else {
                F_right = cell.flux_internal[j];
            }
            
            out << setprecision(10)
                << i << ","
                << j << ","
                << sc.center_x << ","
                << sc.convar[0] << ","
                << sc.convar_old[0] << ","
                << F_left[0] << "," << F_left[1] << "," << F_left[2] << ","
                << F_right[0] << "," << F_right[1] << "," << F_right[2] << ","
                << cell.flux_high_right[0] << "," << cell.flux_high_right[1] << "," << cell.flux_high_right[2] << ","
                << cell.flux_low_right[0] << "," << cell.flux_low_right[1] << "," << cell.flux_low_right[2] << ","
                << cell.alpha << endl;
        }
    }
    
    out.close();
}

void PrintSubCell1dInfo(const FluidSubCell1d& cell, bool verbose) {
    cout << "========================================" << endl;
    cout << "Macro cell index: " << cell.macro_idx << endl;
    cout << "Alpha: " << cell.alpha << " (raw: " << cell.alpha_raw << ")" << endl;
    
    if (verbose) {
        cout << "Legendre coefficients: ";
        for (int j = 0; j <= 4; j++) {
            cout << cell.legendre_coeff[j] << " ";
        }
        cout << endl;
        
        cout << "Sub-cells:" << endl;
        for (int j = 0; j < N_SUBCELL; j++) {
            const SubCell1d& sc = cell.sub[j];
            cout << "  [" << j << "] x=" << sc.center_x 
                 << " dx=" << sc.dx
                 << " w=" << sc.weight
                 << " rho=" << sc.convar[0] << endl;
        }
        
        cout << "Macro fluxes:" << endl;
        cout << "  Left:  " << cell.flux_left[0] << " " << cell.flux_left[1] << " " << cell.flux_left[2] << endl;
        cout << "  Right: " << cell.flux_right[0] << " " << cell.flux_right[1] << " " << cell.flux_right[2] << endl;
    }
    cout << "========================================" << endl;
}
