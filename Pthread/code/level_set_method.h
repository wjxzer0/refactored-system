#ifndef _LEVEL_SET_METHOD_H_
#define _LEVEL_SET_METHOD_H_

#include "Global.hpp"
#include "Matrix.hpp"
#include "Variable.h"

class Clevel_set_method {
public:
  Clevel_set_method();

  ~Clevel_set_method();

  Clevel_set_method(double left, double right, double bottom, double top,
                    int xNum, int yNum, double time, double time_step,
                    CMatrix<CVariable> &Conserved,
                    CMatrix<double> &distance_phi);

public:
  //初始化
  void Initial();

  //空间求解器
  void SpatSolver();

  //返回x方向的距离函数梯度
  CMatrix<double>&  X_dericition_delta_phi();

  //返回y方向的距离函数梯度
  CMatrix<double>&  Y_dericition_delta_phi();

  //设置边界条件
  void level_set_SetBoundary();

  // 5阶J-S WENO
  static double phi_weno_5th(double v1, double v2, double v3, double v4,
                             double v5);

  // x方向（i,j）处速度为正时，重构距离函数梯度
  void X_Weno_plus(int i, int j);

  // x方向（i,j）处速度为负时，重构距离函数梯度
  void X_Weno_mins(int i, int j);

  // y方向（i,j）处速度为正时，重构距离函数梯度
  void Y_Weno_plus(int i, int j);

  // y方向（i,j）处速度为负时，重构距离函数梯度
  void Y_Weno_mins(int i, int j);

  //重新初始化预处理
  double pre_reinitialization(int i, int j);

  //每更新完距离函数就进行一次初始化
  void reinitialization();

  //时间求解器
  void TempSolver(int timepre);

  //输出函数
 CMatrix<double> level_set_Output();

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

  //时间步长
  double m_time_step;

  //速度场
  CMatrix<double> ini_x_velocity;
  CMatrix<double> ini_y_velocity;

  // x方向总网格数+左边界网格数
  int xEnd;

  // x方向总网格数+下边界网格数
  int yEnd;

  // x方向总网格数+左右边界网格数
  int xTotal;

  // y方向总网格数+上下边界网格数
  int yTotal;

  //二维网格
  CMatrix<Center> Grid;

  //距离函数
  CMatrix<double> phi;
  CMatrix<double> phi_temp;
  CMatrix<double> ini_phi;

  // x方向距离函数的导数
  CMatrix<double> X_partial_phi;

  // y方向距离函数的导数
  CMatrix<double> Y_partial_phi;

  //速度场
  CMatrix<Center> velocity;
};
#endif