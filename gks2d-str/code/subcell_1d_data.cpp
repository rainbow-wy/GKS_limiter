/**
 * @file subcell_1d_data.cpp
 * @brief 1D子单元数据结构实现
 * 
 * @note 遵守Zero-Modification原则
 */

#include "subcell_1d_data.h"
#include <cstring>
#include <cmath>

//=============================================================================
// 初始化函数
//=============================================================================

FluidSubCell1d* AllocateSubCell1dArray(Block1d& block) {
    int total = block.nx;
    FluidSubCell1d* fluids_sc = new FluidSubCell1d[total];
    std::memset(fluids_sc, 0, total * sizeof(FluidSubCell1d));
    return fluids_sc;
}

void FreeSubCell1dArray(FluidSubCell1d* fluids_sc) {
    if (fluids_sc != nullptr) {
        delete[] fluids_sc;
    }
}

void InitSubCell1dGeometry(
    FluidSubCell1d& cell,
    Fluid1d* parent,
    int macro_idx
) {
    cell.parent = parent;
    cell.macro_idx = macro_idx;
    
    double macro_dx = parent->dx;
    double macro_cx = parent->cx;
    
    // 初始化光滑度因子
    cell.alpha = 0.0;
    cell.alpha_raw = 0.0;
    for (int j = 0; j <= N_SUBCELL; j++) {
        cell.legendre_coeff[j] = 0.0;
    }
    
    // 初始化通量
    Zero3(cell.flux_left);
    Zero3(cell.flux_right);
    Zero3(cell.flux_high_left);
    Zero3(cell.flux_low_left);
    Zero3(cell.flux_high_right);
    Zero3(cell.flux_low_right);
    for (int i = 0; i < N_SUBCELL - 1; i++) {
        Zero3(cell.flux_internal[i]);
    }
    
    //=========================================================================
    // 核心：基于GL权重的非均匀子单元划分
    // sub_dx = macro_dx * weight[i] / 2.0
    // (除以2是因为参考区间 ξ∈[-1,1] 长度为2)
    //=========================================================================
    
    // 计算子单元左边界（宏观单元左边界）
    double x_left = macro_cx - 0.5 * macro_dx;
    
    for (int i = 0; i < N_SUBCELL; i++) {
        SubCell1d& sc = cell.sub[i];
        
        // 几何属性
        sc.weight = GL_WEIGHTS_5[i];
        sc.dx = macro_dx * sc.weight / 2.0;  // 关键公式！
        
        // 子单元中心位置
        sc.center_x = x_left + 0.5 * sc.dx;
        x_left += sc.dx;  // 移动到下一个子单元左边界
        
        sc.local_idx = i;
        
        // 从宏观单元复制初始值
        Copy3(sc.convar, parent->convar);
        Copy3(sc.convar_old, parent->convar);
        Zero3(sc.primvar);
    }
}

void InitAllSubCells1d(
    FluidSubCell1d* fluids_sc,
    Fluid1d* fluids,
    Block1d& block
) {
    for (int i = 0; i < block.nx; i++) {
        InitSubCell1dGeometry(fluids_sc[i], &fluids[i], i);
    }
}

//=============================================================================
// 数据操作函数
//=============================================================================

void DistributeToSubCells1d(FluidSubCell1d& cell) {
    for (int i = 0; i < N_SUBCELL; i++) {
        Copy3(cell.sub[i].convar, cell.parent->convar);
        Copy3(cell.sub[i].convar_old, cell.parent->convar);
    }
    
}

void ProjectToMacro1d(FluidSubCell1d& cell) {
    // GL加权平均: W_macro = Σ (weight[i]/2) * W_sub[i]
    // 注意：GL权重之和为2，所以除以2得到正确的平均
    double w_macro[3] = {0.0, 0.0, 0.0};
    
    for (int i = 0; i < N_SUBCELL; i++) {
        double w = cell.sub[i].weight / 2.0;
        for (int v = 0; v < 3; v++) {
            w_macro[v] += w * cell.sub[i].convar[v];
        }
    }
    
    Copy3(cell.parent->convar, w_macro);
    Convar_to_primvar_1D(cell.parent->primvar, cell.parent->convar); // 同步更新原始变量
}

void CopySubCell1dNewToOld(FluidSubCell1d& cell) {
    for (int i = 0; i < N_SUBCELL; i++) {
        Copy3(cell.sub[i].convar_old, cell.sub[i].convar);
    }
}

void CopyAllSubCells1dNewToOld(FluidSubCell1d* fluids_sc, Block1d& block) {
    for (int i = 0; i < block.nx; i++) {
        CopySubCell1dNewToOld(fluids_sc[i]);
    }
}

void SubCell1dConvarToPrimvar(FluidSubCell1d& cell) {
    for (int i = 0; i < N_SUBCELL; i++) {
        SubCell1d& sc = cell.sub[i];
        double rho = sc.convar[0];
        
        if (rho <= 0.0) {
            sc.primvar[0] = rho;
            sc.primvar[1] = 0.0;
            sc.primvar[2] = 0.0;
            continue;
        }
        
        double u = sc.convar[1] / rho;
        double E = sc.convar[2];
        double p = (Gamma - 1.0) * (E - 0.5 * rho * u * u);
        
        sc.primvar[0] = rho;
        sc.primvar[1] = u;
        sc.primvar[2] = p;
    }
}

bool IsAdmissible1d(const double* convar, double eps) {
    double rho = convar[0];
    if (rho <= eps || std::isnan(rho)) {
        return false;
    }
    
    double rho_u = convar[1];
    double E = convar[2];
    double u = rho_u / rho;
    double p = (Gamma - 1.0) * (E - 0.5 * rho * u * u);
    
    if (p <= eps || std::isnan(p)) {
        return false;
    }
    
    return true;
}
