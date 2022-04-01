#ifndef _VARIABLE_H_
#define _VARIABLE_H_
#pragma once

#include <cassert>
#include <cmath>

using namespace std;

class CVariable {
public:
  /**
   * @brief 默认构造函数
   *
   */
  CVariable();

  /**
   * @brief 构造函数，构造4x1的列向量
   *
   * @param a double型变量
   * @param b double型变量
   * @param c double型变量
   * @param d double型变量
   */
  CVariable(double a, double b, double c, double d);

  //~CVariable();
public:
  /**
   * @brief 向量的四个分量
   *
   */
  double u1;
  double u2;
  double u3;
  double u4;

public:
  /**
   * @brief 向量赋值
   *
   * @param v CVariable型变量
   * CVariable& 引用
   */
  CVariable &operator=(const CVariable &v);

  /**
   * @brief 向量加法
   *
   * @param v  CVariable型变量
   * CVariable 返回结果
   */
  CVariable operator+(const CVariable &v);

  /**
   * @brief 向量减法
   *
   * @param v CVariable型变量
   * CVariable 返回结果
   */
  CVariable operator-(const CVariable &v);

  /**
   * @brief 向量左数乘
   *
   * @param c double型变量
   * @param v CVariable型变量
   * CVariable 返回结果
   */
  friend CVariable operator*(double c, const CVariable &v);

  /**
   * @brief 向量内积
   *
   * @param v1 CVariable型变量
   * @param v2 CVariable型变量
   * double 返回结果
   */
  friend double operator*(const CVariable &v1, const CVariable &v2);

  /**
   * @brief 向量中取最大绝对值
   *
   * double 返回结果
   */
  double GetAbsMax();
};




//下面的变量是为了加WENO高阶重构定义的，只保留着未使用
class CPhyVar {
public:
  CPhyVar();

public:
  double rho;
  double u;
  double v;
  double p;
  double e;
  double H;
  double a; // sound speed
};

/**
 * @brief 定义一个数据结构，用于储存网格信息
 *
 */
struct Center {
  double x;
  double y;
};

/**
 * @brief 定义一个数据结构，用于储存计算区域信息
 *
 */
struct Rect {
  double Left;
  double Right;
  double Top;
  double Bottom;
};
/**
 * @brief 定义一个数据结构，用于储存WENO的非线性权
 *
 */
struct Tri {
  double beta1;
  double beta2;
  double beta3;
};
#endif