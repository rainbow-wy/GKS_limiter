/**
 * @file riemann_problem_subcell_1d.cpp
 * @brief 1D子单元GKS求解器测试 - 黎曼问题配置
 * 
 * =============================================================================
 *                         算例配置说明
 * =============================================================================
 * 
 * 1. 在下方"算例选择"区域，取消注释您想运行的算例
 * 2. 确保只有一个算例未被注释
 * 3. 重新编译并运行
 * 
 * =============================================================================
 */

#include "riemann_problem_subcell_1d.h"
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

//=============================================================================
//                              算例选择
// 
// ★★★ 取消注释您想运行的算例 ★★★
//=============================================================================

//-----------------------------------------------------------------------------
// 【算例 1: Double Rarefaction (双稀疏波)】
// 描述：测试近真空状态下的稳定性，包含两个向外扩散的稀疏波
// 推荐参数：CFL=0.5, Gamma=1.4
//-----------------------------------------------------------------------------
//#define CASE_DOUBLE_RAREFACTION

//-----------------------------------------------------------------------------
// 【算例 2: Leblanc Shock Tube (Leblanc 激波管)】
// 描述：极端激波测试，压力比高达 10^9
// 推荐参数：CFL=0.3 (需要更小的时间步), Gamma=1.4
// 注意：此算例需要极高的时间步长精度，建议使用更小的CFL数
//-----------------------------------------------------------------------------
//#define CASE_LEBLANC

//-----------------------------------------------------------------------------
// 【算例 3: Sod Shock Tube (Sod 激波管)】 - 经典算例
// 描述：经典Sod问题，测试激波、接触间断和稀疏波
// 推荐参数：CFL=0.5, Gamma=1.4
//-----------------------------------------------------------------------------
#define CASE_SOD

//-----------------------------------------------------------------------------
// 【算例 4: Lax Shock Tube (Lax 激波管)】
// 描述：Lax问题，比Sod更强的激波
// 推荐参数：CFL=0.5, Gamma=1.4
//-----------------------------------------------------------------------------
//#define CASE_LAX

//-----------------------------------------------------------------------------
// 【算例 5: Shu-Osher (激波与正弦波相互作用)】
// 描述：激波与密度扰动的相互作用
// 推荐参数：CFL=0.5, Gamma=1.4
//-----------------------------------------------------------------------------
//#define CASE_SHU_OSHER

//=============================================================================
//                           算例参数定义
//=============================================================================

#ifdef CASE_DOUBLE_RAREFACTION
//-----------------------------------------------------------------------------
// 双稀疏波 (Double Rarefaction)
//-----------------------------------------------------------------------------
const double X_LEFT    = -1.0;      // 左边界
const double X_RIGHT   = 1.0;       // 右边界
const int    MESH_NUM  = 200;       // 网格数
const double T_END     = 0.6;       // 结束时间
const double CFL_NUM   = 0.5;       // CFL数
const double GAMMA_GAS = 1.4;       // 比热比
const double X_DISCON  = 0.0;       // 间断位置

// 初始条件 (rho, u, p)
const double RHO_L = 7.0,  U_L = -1.0, P_L = 0.2;   // 左侧状态
const double RHO_R = 7.0,  U_R =  1.0, P_R = 0.2;   // 右侧状态

const char* CASE_NAME = "Double Rarefaction";
#endif

#ifdef CASE_LEBLANC
//-----------------------------------------------------------------------------
// Leblanc 激波管 (极端压力比 10^9)
// 警告：此算例对数值格式要求极高！
//-----------------------------------------------------------------------------
const double X_LEFT    = -1.0;      // 左边界
const double X_RIGHT   = 1.0;       // 右边界
const int    MESH_NUM  = 200;       // 网格数
const double T_END     = 0.0001;    // 结束时间 (很短，因为压力差极大)
const double CFL_NUM   = 0.3;       // CFL数 (需要更小)
const double GAMMA_GAS = 1.4;       // 比热比
const double X_DISCON  = 0.0;       // 间断位置

// 初始条件 (rho, u, p)
const double RHO_L = 2.0,   U_L = 0.0, P_L = 1.0e9;  // 左侧状态 (高压)
const double RHO_R = 0.001, U_R = 0.0, P_R = 1.0;    // 右侧状态 (低压)

const char* CASE_NAME = "Leblanc Shock Tube (Extreme)";
#endif

#ifdef CASE_SOD
//-----------------------------------------------------------------------------
// 经典 Sod 激波管
//-----------------------------------------------------------------------------
const double X_LEFT    = 0.0;       // 左边界
const double X_RIGHT   = 1.0;       // 右边界
const int    MESH_NUM  = 200;       // 网格数
const double T_END     = 0.2;       // 结束时间
const double CFL_NUM   = 0.05;       // CFL数
const double GAMMA_GAS = 1.4;       // 比热比
const double X_DISCON  = 0.5;       // 间断位置

// 初始条件 (rho, u, p)
const double RHO_L = 1.0,   U_L = 0.0, P_L = 1.0;   // 左侧状态
const double RHO_R = 0.125, U_R = 0.0, P_R = 0.1;   // 右侧状态

const char* CASE_NAME = "Sod Shock Tube";
#endif

#ifdef CASE_LAX
//-----------------------------------------------------------------------------
// Lax 激波管
//-----------------------------------------------------------------------------
const double X_LEFT    = 0.0;       // 左边界
const double X_RIGHT   = 1.0;       // 右边界
const int    MESH_NUM  = 200;       // 网格数
const double T_END     = 0.15;      // 结束时间
const double CFL_NUM   = 0.5;       // CFL数
const double GAMMA_GAS = 1.4;       // 比热比
const double X_DISCON  = 0.5;       // 间断位置

// 初始条件 (rho, u, p)
const double RHO_L = 0.445, U_L = 0.698, P_L = 3.528;  // 左侧状态
const double RHO_R = 0.5,   U_R = 0.0,   P_R = 0.571;  // 右侧状态

const char* CASE_NAME = "Lax Shock Tube";
#endif

#ifdef CASE_SHU_OSHER
//-----------------------------------------------------------------------------
// Shu-Osher 问题 (激波与正弦波相互作用)
// 注意：这个算例的初始化需要特殊处理（包含正弦扰动）
//-----------------------------------------------------------------------------
const double X_LEFT    = -5.0;      // 左边界
const double X_RIGHT   = 5.0;       // 右边界
const int    MESH_NUM  = 400;       // 网格数 (需要更多)
const double T_END     = 1.8;       // 结束时间
const double CFL_NUM   = 0.5;       // CFL数
const double GAMMA_GAS = 1.4;       // 比热比
const double X_DISCON  = -4.0;      // 激波初始位置

// 初始条件 (rho, u, p) - 左侧为激波后状态
const double RHO_L = 3.857143, U_L = 2.629369, P_L = 10.33333;  // 激波后
const double RHO_R = 1.0,      U_R = 0.0,      P_R = 1.0;       // 激波前 (基值)

const char* CASE_NAME = "Shu-Osher Problem";
#endif

//=============================================================================
//                           初始条件函数
//=============================================================================

void ICfor1dSubcell(
    Fluid1d* fluids,
    double rho_L, double u_L, double p_L,
    double rho_R, double u_R, double p_R,
    double x_discon,
    Block1d& block
) {
    for (int i = 0; i < block.nx; i++) {
        double x = fluids[i].cx;
        double rho, u, p;
        
#ifdef CASE_SHU_OSHER
        // Shu-Osher 特殊初始化：激波前有正弦扰动
        if (x < x_discon) {
            rho = rho_L;
            u = u_L;
            p = p_L;
        } else {
            // 右侧密度有正弦扰动
            rho = rho_R + 0.2 * sin(5.0 * x);
            u = u_R;
            p = p_R;
        }
#else
        // 标准 Riemann 问题初始化
        if (x <= x_discon) {
            rho = rho_L;
            u = u_L;
            p = p_L;
        } else {
            rho = rho_R;
            u = u_R;
            p = p_R;
        }
#endif
        
        // 原始变量
        fluids[i].primvar[0] = rho;
        fluids[i].primvar[1] = u;
        fluids[i].primvar[2] = p;
        
        // 守恒变量
        fluids[i].convar[0] = rho;
        fluids[i].convar[1] = rho * u;
        fluids[i].convar[2] = p / (Gamma - 1.0) + 0.5 * rho * u * u;
    }
}

//=============================================================================
//                           主测试函数
//=============================================================================

void riemann_problem_subcell_1d() {
    cout << "============================================" << endl;
    cout << "  1D Sub-cell GKS Solver" << endl;
    cout << "  Test Case: " << CASE_NAME << endl;
    cout << "============================================" << endl;
    
    // 计时
    clock_t time_start = clock();
    
    //=========================================================================
    // 设置全局参数
    //=========================================================================
    Gamma = GAMMA_GAS;
    K = 4;  // 1D: K = (5 - 3*Gamma) / (Gamma - 1) 对于单原子气体
    R_gas = 1.0;
    tau_type = Euler;
    c1_euler = 0.05;
    c2_euler = 1.0;
    
    // 设置默认重构方法
    cellreconstruction = Vanleer;
    wenotype = wenoz;
    reconstruction_variable = characteristic;
    g0reconstruction = Center_3rd;
    
    //=========================================================================
    // 初始化 Block
    //=========================================================================
    Block1d block;
    block.ghost = 3;           // WENO5需要3个ghost cell
    block.nodex = MESH_NUM;
    block.nx = MESH_NUM + 2 * block.ghost;
    block.left = X_LEFT;
    block.right = X_RIGHT;
    block.dx = (X_RIGHT - X_LEFT) / MESH_NUM;
    block.CFL = CFL_NUM;
    block.t = 0.0;
    block.step = 0;
    block.stages = 2;
    
    // 初始化时间推进系数
    Initial_stages(block);
    
    cout << "\nMesh Configuration:" << endl;
    cout << "  Domain: [" << X_LEFT << ", " << X_RIGHT << "]" << endl;
    cout << "  Cells: " << MESH_NUM << endl;
    cout << "  dx: " << block.dx << endl;
    cout << "  CFL: " << CFL_NUM << endl;
    cout << "  T_end: " << T_END << endl;
    cout << "  Gamma: " << GAMMA_GAS << endl;
    
    //=========================================================================
    // 分配宏观单元数组
    //=========================================================================
    Fluid1d* fluids = new Fluid1d[block.nx];
    
    // 设置网格坐标
    for (int i = 0; i < block.nx; i++) {
        fluids[i].dx = block.dx;
        fluids[i].cx = block.left + (i + 0.5 - block.ghost) * block.dx;
    }
    
    //=========================================================================
    // 设置初始条件
    //=========================================================================
    cout << "\nInitial Condition:" << endl;
    cout << "  Left:  (rho=" << RHO_L << ", u=" << U_L << ", p=" << P_L << ")" << endl;
    cout << "  Right: (rho=" << RHO_R << ", u=" << U_R << ", p=" << P_R << ")" << endl;
    cout << "  Discontinuity at x = " << X_DISCON << endl;
    
    ICfor1dSubcell(fluids, RHO_L, U_L, P_L, RHO_R, U_R, P_R, X_DISCON, block);
    
    //=========================================================================
    // 边界条件函数指针
    //=========================================================================
    BoundaryCondition leftboundary = free_boundary_left;
    BoundaryCondition rightboundary = free_boundary_right;
    Fluid1d bcvalue[2];  // 边界值 (透射边界不需要具体值)
    
    //=========================================================================
    // 初始化子单元求解器
    //=========================================================================
    FluidSubCell1d* fluids_sc = InitSubCell1dSolver(block, fluids);
    
    // 初始输出
    OutputSubCell1d(fluids_sc, block);
    
    //=========================================================================
    // 时间推进
    //=========================================================================
    cout << "\n" << "Starting time evolution..." << endl;
    cout << "----------------------------------------" << endl;
    
    clock_t time_compute_start = clock();
    
    int output_interval = 100;  // 输出间隔
    
    while (block.t < T_END) {
        // 保存旧状态
        CopyAllSubCells1dNewToOld(fluids_sc, block);
        
        // 计算时间步长
        block.dt = Get_CFL(block, fluids, T_END);
        
        // 边界条件 (透射边界)
        leftboundary(fluids, block, bcvalue[0]);
        rightboundary(fluids, block, bcvalue[1]);
        
        // 传递边界到子单元
        ApplySubCell1dBC(fluids_sc, fluids, block);
        
        // 子单元求解器推进一步
        SubCell1dSolverStep(fluids_sc, fluids, block);
        
        // 更新时间
        block.t += block.dt;
        block.step++;
        
        // 定期输出 (已根据用户要求禁用中间输出)
        /*
        if (block.step % output_interval == 0) {
            OutputSubCell1d(fluids_sc, block);
        }
        */
        
        // 第一步后输出通量调试信息
        if (block.step == 1) {
            OutputFluxDebug(fluids_sc, block);
        }
        
        // 检查是否爆炸
        if (fluids[block.nx / 2].convar[0] != fluids[block.nx / 2].convar[0]) {
            cout << "ERROR: Simulation blows up!" << endl;
            OutputSubCell1d(fluids_sc, block);
            break;
        }
    }
    
    //=========================================================================
    // 最终输出
    //=========================================================================
    OutputSubCell1d(fluids_sc, block);
    
    clock_t time_end = clock();
    
    cout << "----------------------------------------" << endl;
    cout << "Simulation Completed!" << endl;
    cout << "  Final time: " << block.t << endl;
    cout << "  Total steps: " << block.step << endl;
    cout << "  Total time: " << (double)(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl;
    cout << "  Compute time: " << (double)(time_end - time_compute_start) / CLOCKS_PER_SEC << " seconds" << endl;
    
    //=========================================================================
    // 清理
    //=========================================================================
    FinalizeSubCell1dSolver(fluids_sc);
    delete[] fluids;
}
