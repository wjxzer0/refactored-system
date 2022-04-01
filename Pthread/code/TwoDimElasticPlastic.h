/********************************************
// File name: TwoDimElasticPlastic.h
//
// Abstract:  This is the C++ source code
//            for Two dimension Euler equations of elastic plastic flow.
//
// Description of some key points:
//
//            1.Only suit for uniform grid;
//            2.The initial data must be constant in the computational area;
//            3.SSHLL/UPWIND/HLLC Riemann solver(HLL Type) are used;
//
//
//
//
// Current version: 1.0
// Author: Zeng zhiqiang(1257371357@qq.com)
// Date: September,21 2020
***********************************************/

#include "Approximate_Riemann_solver.hpp"
#include "Global.hpp"
#include "Matrix.hpp"
//#include "MyMath.hpp"
#include "Variable.h"
#include "level_set_function.hpp"
#include"level_set_method.h"
#include"cell.hpp"


using namespace std;

class CTwoDimElasticPlastic {
public:
  /**
   * @brief 默认构造函数
   *
   */
  CTwoDimElasticPlastic();

  /**
   * @brief 析构函数
   *
   */
  ~CTwoDimElasticPlastic();

  /**
   * @brief 构造弹塑性流场计算
   *
   * @param left 左边界
   * @param right 右边界
   * @param bottom 下边界
   * @param top 上边界
   * @param xNum x方向网格数
   * @param yNum y方向网格数
   * @param time 总时间
   * @param leftstate 左初值，包含rho , u ,v, p
   * @param rightstate 右初值，包含rho , u ,v, p
   * @param leftstress 左初值，包含s_xx,s_yy,s_xy
   * @param rightstress 右初值，包含s_xx,s_yy,s_xy
   * @param y_e 初始塑性屈服强度
   */
  CTwoDimElasticPlastic(double left, double right, double bottom, double top,
                        int xNum, int yNum, double time,
                        vector<double> &leftstate, vector<double> &rightstate,
                        vector<double> &leftstress, vector<double> &rightstress,
                        double y_e);

public:
  /**
   * @brief 平面波算例初始化函数
   */
  void Initial1();

  /**
   * @brief 钉钉子算例初始化函数
   */
  void Initial2();

  /**
   * @brief 边界赋值,定义了三个虚拟网格
   *
   * @param time 时间
   */
  void SetBoundary(double time);

public:
  /**
   * @brief 获取x方向数值通量
   *
   */
  void GetX_flux();

  /**
   * @brief 获取x方向星区域速度
   *
   */
  void GetX_Uspeed();

public:
  /**
   * @brief 获取y方向数值通量
   *
   */
  void GetY_flux();

  /**
   * @brief 获取y方向星区域速度
   *
   */
  void GetY_Uspeed();

public:
  /**
   * @brief 获取时间步长
   *
   * double 返回时间步长
   */
  double GetTimeStep(void);

  /**
   * @brief 时间推进
   *
   * @param t n时刻的时间
   */
  void Timestep(double t);

  /**
   * @brief 3rd Runge-kutta方法，守恒物理量三阶，应力一阶
   *
   * @param t n时刻的时间
   * @param m_tstep 时间步长
   */
  void Rungekutta(double t, double m_tstep);

  /**
   * @brief 主计算函数
   *
   */
  void MainSolve();

  /**
   * @brief 输出
   *
   */
  void Output();

  // void Real_area_index(double time, int x_point, int point_y);

  void Virtual_assignment(double time);

private:
  //计算区域左值
  double m_left;

  //计算区域右值
  double m_right;

  //计算区域上值
  double m_top;

  //计算区域下值
  double m_bottom;

  // x方向总网格数
  int m_xNum;

  // y方向总网格数
  int m_yNum;

  // x方向空间步长
  double m_xH;

  // y方向空间步长
  double m_yH;

  //计算终止时间
  double m_Time;

  // x方向总网格数+左边界网格数
  int xEnd;

  // y方向总网格数+下边界网格数
  int yEnd;

  // x方向总网格数+左右边界网格数
  int xTotal;

  // y方向总网格数+上下边界网格数
  int yTotal;

  //初始材料屈服强度
  double m_ini_Y;

  //初始左密度
  double m_ini_lrho;

  //初始左速度u
  double m_ini_lu;

  //初始左速度v
  double m_ini_lv;

  //初始左压强
  double m_ini_lp;

  //初始左偏应力s_xx
  double m_ini_lsx;

  //初始左偏应力s_yy
  double m_ini_lsy;

  //初始左偏应力s_xy
  double m_ini_lsxy;

  //初始右密度
  double m_ini_rrho;

  //初始右速度u
  double m_ini_ru;

  //初始右速度v
  double m_ini_rv;

  //初始右压强
  double m_ini_rp;

  //初始右偏应力s_xx
  double m_ini_rsx;

  //初始右偏应力s_yy
  double m_ini_rsy;

  //初始右偏应力s_xy
  double m_ini_rsxy;

  //二维网格
  CMatrix<Center> Grid;

  //守恒物理量
  CMatrix<CVariable> U0;

  // x方向守恒通量
  CMatrix<CVariable> FU0;

  // y方向守恒通量
  CMatrix<CVariable> GU0;

  //应力向量
  CMatrix<CVariable> S0;

  //线性硬化屈服应力
  CMatrix<CVariable> Y_E;

  //线性硬化判断系数
  CMatrix<CVariable> tempY;

  // x方向一维通量
  CMatrix<CVariable> Fflux;

  // y方向一维通量
  CMatrix<CVariable> Gflux;

  // y方向星区域速度v
  CMatrix<CVariable> Y_speed;

  // y方向星区域速度u
  CMatrix<CVariable> u_Y_speed;

  // x方向星区域的速度u
  CMatrix<CVariable> X_speed;

  // x方向的星区域速度v
  CMatrix<CVariable> v_X_speed;

  //时间推进中守恒量的中间变量
  CMatrix<CVariable> tempU0;
  CMatrix<CVariable> tempU1;
  CMatrix<CVariable> tempU2;
  CMatrix<CVariable> tempU3;

  //时间推进中应力的中间变量
  CMatrix<CVariable> tempS0;
  CMatrix<CVariable> tempS1;
  CMatrix<CVariable> tempS2;
  CMatrix<CVariable> tempS3;

  //时间推进中线性硬化屈服强度的中间变量
  CMatrix<CVariable> tempY0;
  CMatrix<CVariable> tempY1;
  CMatrix<CVariable> tempY2;
  CMatrix<CVariable> tempY3;

   //距离函数
    CMatrix<double>  distance_phi;
    CMatrix<double> distance_phi_temp;

    // //界面外虚拟赋值的守恒物理量和应力向量
    // CMatrix<CVariable> left_U0;
    // CMatrix<CVariable> right_U0;
    // CMatrix<CVariable> top_U0;
    // CMatrix<CVariable> bottom_U0;
  
    // CMatrix<CVariable> left_S0;
    // CMatrix<CVariable> right_S0;
    // CMatrix<CVariable> top_S0;
    // CMatrix<CVariable> bottom_S0;


     // x方向距离函数的导数
  CMatrix<double> X_dericition_delta_phi;

  // y方向距离函数的导数
  CMatrix<double> Y_dericition_delta_phi;


};