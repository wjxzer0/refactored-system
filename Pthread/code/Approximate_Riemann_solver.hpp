/******************************************
 * Filename: Approximate_Riemann_solver.hpp
 * File title: Two dimensional elastic-plastic approximate Riemann solver
 * Abstract: This is a function library containing upwind, HLL and HLLC
 *approximate Riemannian solvers, which can be directly used to solve
 *two-dimensional elastic-plastic Riemannian problems
 *
 * Modification date:	January 23, 2020
 * Author:ZZQ_BUAA
 ********************************************/

#ifndef _2DEulerFun_H_
#define _2DEulerFun_H_
#pragma once

#include "Global.hpp"
#include "Variable.h"

using namespace std;

/**
 * @brief 根据状态方程获取内能
 *
 * @param rho 密度
 * @param p 压强
 * double 返回内能
 */
inline double internal_energy(const double rho, const double p) {
  return (p / Gamma1 - const_c * const_c * (rho - const_rho) / Gamma1);
}

/**
 * @brief 获得总能量
 *
 * @param rho 密度
 * @param u  速度
 * @param v  速度
 * @param p  压强
 * double 返回总能量
 */
inline double energy(const double rho, const double u, const double v,
                     const double p) {
  return (0.5 * rho * u * u + 0.5 * rho * v * v + internal_energy(rho, p));
}

/**
 * @brief 根据状态方程获取压力p
 *
 * @param U 守恒向量
 * double 返回压力p
 */
inline double obtain_p(CVariable &U) {
  return (Gamma1 * U.u4 - 0.5 * Gamma1 * (U.u2 * U.u2 / U.u1) -
          0.5 * Gamma1 * (U.u3 * U.u3 / U.u1) +
          const_c * const_c * (U.u1 - const_rho));
}

/**
 * @brief 根据状态方程获取声速
 *
 * @param U 守恒物理量
 * @param s 偏应力
 * double 返回声速
 */
inline double sonic_speed(CVariable &U, double s) {
  double rho(0), p(0), u(0);

  rho = U.u1;
  u = U.u2 / U.u1;
  p = obtain_p(U);

  return (
      sqrt((Gamma * p + const_c * const_c * const_rho - Gamma1 * s) / (rho)));
}

/**
 * @brief 获得x方向的通量
 *
 * @param u 守恒量
 * @param s
 * 应力，第一分量为x方向偏应力s_xx，第二分量为y方向偏应力s_yy，第三分量为剪切力s_xy，第四分量恒为0
 * CVariable 返回x方向通量
 */
inline CVariable FU(CVariable &u, const CVariable &s) {
  double s0(0), s1(0), s2(0), s3(0), vy(0), p(0);

  p = obtain_p(u);
  vy = u.u3 / u.u1;

  s0 = u.u2;
  s1 = u.u2 * u.u2 / u.u1 + p - s.u1; // x方向通量的四个分量
  s2 = u.u2 * u.u3 / u.u1 - s.u3;
  s3 = u.u2 / u.u1 * (u.u4 + p - s.u1) - vy * s.u3;

  return CVariable(s0, s1, s2, s3);
}

/**
 * @brief 获得y方向的通量
 *
 * @param u 守恒量
 * @param s
 * 应力，第一分量为x方向偏应力s_xx，第二分量为y方向偏应力s_yy，第三分量为剪切力s_xy，第四分量恒为0
 * CVariable 返回y方向通量
 */
inline CVariable GU(CVariable &u, const CVariable &s) {
  double s0(0), s1(0), s2(0), s3(0), ux(0), p(0);

  p = obtain_p(u);
  ux = u.u2 / u.u1;

  s0 = u.u3;
  s1 = u.u2 * u.u3 / u.u1 - s.u3; // y方向通量的四个分量
  s2 = u.u3 * u.u3 / u.u1 + p - s.u2;
  s3 = u.u3 / u.u1 * (u.u4 + p - s.u2) - ux * s.u3;

  return CVariable(s0, s1, s2, s3);
}

/**
 * @brief 判断偏应力是否满足米塞斯屈服极限并拉回屈服极限
 *
 * @param a
 * 应力，第一分量为x方向偏应力s_xx，第二分量为y方向偏应力s_yy，第三分量为剪切力s_xy，第四分量恒为0
 * CVariable 返回修正后的应力
 */
inline CVariable Limit(CVariable &a) {
  CVariable mid;
  double M;
  M = (pow(abs(a.u1), 2.0) + pow(abs(a.u2), 2.0) + 2.0 * pow(abs(a.u3), 2.0) +
       pow(abs(a.u1 + a.u2), 2.0)) /
      ((2.0 / 3.0) * pow(Y0, 2.0)); //判断是否超过米塞斯屈服极限
  if (abs(M - 1) > Eps) {
    mid = sqrt(1.0 / M) * a;
  } else {
    mid = a;
  }
  return (mid);
}

/**
 * @brief 判断是否超过米塞斯屈服极限
 *
 * @param a
 * 应力，第一分量为x方向偏应力s_xx，第二分量为y方向偏应力s_yy，第三分量为剪切力s_xy，第四分量恒为0
 * double 返回屈服系数
 */
inline double vos_limit(CVariable &a) {
  CVariable mid;
  double M;

  M = (pow(abs(a.u1), 2.0) + pow(abs(a.u2), 2.0) + 2.0 * pow(abs(a.u3), 2.0) +
       pow(abs(a.u1 + a.u2), 2.0)) /
      ((2.0 / 3.0) * pow(Y0, 2.0)); //给出应力与屈服极限之间的比例

  return (M);
}

/**
 * @brief HLL近似解法器计算x方向星区域的偏应力
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回星区域应力，用于判断是否超过米塞斯屈服极限
 */
inline CVariable stress_x(CVariable &a, CVariable &b, CVariable &s1,
                          CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0);

  CVariable U_S, stress_s, mid;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1; //给出原始物理量
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1));
  sl = -sr;

  U_S = 1.0 / (sr - sl) *
        (sr * b - sl * a + FU(a, s1) - FU(b, s2)); //获取星区域守恒物理量

  stress_s.u1 = min((4.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u1,
                    (4.0 / 3.0) * Miu * log(b.u1 / U_S.u1) + s2.u1);
  stress_s.u2 = min(-(2.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u2,
                    -(2.0 / 3.0) * Miu * log(b.u1 / U_S.u1) +
                        s2.u2); //获取星区域应力s_xx,s_yy,s_xy

  stress_s.u3 = min(s1.u3, s2.u3);
  stress_s.u4 = 0.0;

  return (stress_s);
}

/**
 * @brief HLLC近似解法器计算x方向左星区域的偏应力
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回左侧星区域应力
 */
inline CVariable HLLC_stress_xL(CVariable &a, CVariable &b, CVariable &s1,
                                CVariable &s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_s;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1; //给出原始物理量
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u1;
  sigmar = -pr + s2.u1;

  cr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1));
  cl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1)); //估计波速

  Hl = rhol * (cl - ul);
  Hr = rhor * (cr - ur);

  speed_u = (-Hl * ul + Hr * ur - sigmal + sigmar) / (Hr - Hl); //获取星区域波速
  rhols = Hl / (cl - speed_u);
  rhors = Hr / (cr - speed_u); //获取星区域密度

  stress_s.u1 = (4.0 / 3.0) * Miu * log(a.u1 / rhols) + s1.u1;
  stress_s.u2 = -(2.0 / 3.0) * Miu * log(a.u1 / rhols) + s1.u2;
  stress_s.u3 = s1.u3;
  stress_s.u4 = 0.0; //获取星区域应力s_xx,s_yy,s_xy

  return (stress_s);
}

/**
 * @brief HLLC近似解法器计算x方向右星区域的偏应力
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向右侧星区域应力
 */
inline CVariable HLLC_stress_xR(CVariable &a, CVariable &b, CVariable &s1,
                                CVariable &s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_s;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1; //给出原始物理量
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u1;
  sigmar = -pr + s2.u1;

  cr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速
  cl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  Hl = rhol * (cl - ul);
  Hr = rhor * (cr - ur);

  speed_u = (-Hl * ul + Hr * ur - sigmal + sigmar) / (Hr - Hl); //获取星区域波速
  rhols = Hl / (cl - speed_u);
  rhors = Hr / (cr - speed_u); //获取星区域密度

  stress_s.u1 = (4.0 / 3.0) * Miu * log(a.u1 / rhors) + s1.u1;
  stress_s.u2 = -(2.0 / 3.0) * Miu * log(a.u1 / rhors) + s1.u2;
  stress_s.u3 = s1.u3;
  stress_s.u4 = 0.0; //获取星区域应力s_xx,s_yy,s_xy

  return (stress_s);
}

/**
 * @brief
 * HLLC近似解法器计算y方向左星区域的偏应力（算法思路与x方向HLLC_stress_xL是一致的）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向左侧星区域应力
 */
inline CVariable HLLC_stress_yL(CVariable &a, CVariable &b, CVariable &s1,
                                CVariable &s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_s;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u2;
  sigmar = -pr + s2.u2;

  cr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  cl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  Hl = rhol * (cl - vl);
  Hr = rhor * (cr - vr);

  speed_v = (-Hl * vl + Hr * vr - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_v);
  rhors = Hr / (cr - speed_v);

  stress_s.u2 = (4.0 / 3.0) * Miu * log(a.u1 / rhols) + s1.u1;
  stress_s.u1 = -(2.0 / 3.0) * Miu * log(a.u1 / rhols) + s1.u2;
  stress_s.u3 = s1.u3;
  stress_s.u4 = 0.0;

  return (stress_s);
}

/**
 * @brief
 * HLLC近似解法器计算y方向右星区域的偏应力（算法思路与HLLC_stress_xR是一致的）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向右侧星区域应力
 */
inline CVariable HLLC_stress_yR(CVariable &a, CVariable &b, CVariable &s1,
                                CVariable &s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_s;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u2;
  sigmar = -pr + s2.u2;

  cr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  cl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  Hl = rhol * (cl - vl);
  Hr = rhor * (cr - vr);

  speed_v = (-Hl * vl + Hr * vr - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_v);
  rhors = Hr / (cr - speed_v);

  stress_s.u2 = (4.0 / 3.0) * Miu * log(a.u1 / rhors) + s1.u1;
  stress_s.u1 = -(2.0 / 3.0) * Miu * log(a.u1 / rhors) + s1.u2;
  stress_s.u3 = s1.u3;
  stress_s.u4 = 0.0;

  return (stress_s);
}

/**
 * @brief HLL近似解法器计算x方向弹性极限状态下的偏应力
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回返回x方向弹性极限状态下的应力：s_xx,s_yy,s_xy
 */
inline CVariable stress_Limit_x(CVariable &a, CVariable &b, CVariable &s1,
                                CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0);

  CVariable U_S, stress_s, mid;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1; //获取原始物理量
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速
  sl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  U_S = 1.0 / (sr - sl) *
        (sr * b - sl * a + FU(a, s1) - FU(b, s2)); //获取星区域守恒量

  stress_s.u1 = (4.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u1;
  stress_s.u2 = -(2.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u2;
  stress_s.u3 = s1.u3;
  stress_s.u4 = 0.0; //获取星区域应力s_xx,s_yy,s_xy

  mid = Limit(stress_s); //假如超过米塞斯屈服极限，就径向拉回

  return (mid);
}

/**
 * @brief HLL近似解法器计算x方向左弹性极限状态下的守恒物理量
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回返回x方向左侧弹性极限状态下的守恒物理量
 */
inline CVariable U_Limit_x_L(CVariable &a, CVariable &b, CVariable &s1,
                             CVariable &s2) {
  double rho(0), speed(0), ps(0), rhos(0), vs(0), us(0);
  double u(0), v(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double p(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), ts(0);
  double rhor(0), ur(0), pr(0);

  CVariable U_S, stress_s, mid;

  rhor = b.u1;
  ur = b.u2 / b.u1;
  pr = obtain_p(b);

  rho = a.u1;
  u = a.u2 / a.u1; //获取原始物理量
  v = a.u3 / a.u1;

  p = obtain_p(a);

  sr = max(u + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1));
  sl = -sr;

  U_S = 1.0 / (sr - sl) *
        (sr * b - sl * a + FU(a, s1) - FU(b, s2)); //获取星区域守恒物理量

  mid.u1 = (4.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u1;
  mid.u2 = -(2.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u2;
  mid.u3 = s1.u3;
  mid.u4 = 0.0; //估计星区域应力s_xx,s_yy,s_xy

  stress_s = Limit(mid); //判断是否超过米塞斯屈服极限，径向拉回

  rhos = rho * exp(-3.0 * (stress_s.u1 - s1.u1) / (4.0 * Miu));
  if (abs(rhos - rho) < Eps_rho) {
    ps = p;
    sigma_s = -ps + stress_s.u1;
    us = u;
    vs = v;
  } else {
    ts = rhos * rho / (rhos - rho);
    ps = (2.0 * Gamma1 * ts * rhos) / (2.0 * ts - Gamma1 * rhos) *
         (const_c * const_c / Gamma1 *
              ((rhos - const_rho) / rhos - (rho - const_rho) / rho) +
          Gamma1 * p / rho);
    sigma_s = -ps + stress_s.u1; //获取弹性极限状态下的原始物理量
    us = u - sqrt((-p + s1.u1 + ps - stress_s.u1) / ts);
    vs = v;
  }

  U_S.u1 = rhos;
  U_S.u2 = rhos * us;
  U_S.u3 = rhos * vs; //用原始物理量构造守恒物理量
  U_S.u4 = energy(rhos, us, vs, ps);

  return (U_S);
}

/**
 * @brief HLLC近似解法器计算x方向左弹性极限状态下的物理量
 *
 * @param a i-1处守恒变量
 * @param s1 i-1处应力
 * @param mid 预估的左侧星区域的应力
 * CVariable 返回返回x方向左侧弹性极限状态下的守恒物理量
 */
inline CVariable HLLC_Limit_x_L(CVariable &a, CVariable &s1, CVariable &mid) {
  double rho(0), speed(0), ps(0), rhos(0), vs(0), us(0);
  double u(0), v(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double p(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), ts(0);
  double rhor(0), ur(0), pr(0);

  CVariable U_S, stress_s;

  rho = a.u1;
  u = a.u2 / a.u1; //获取原始物理量
  v = a.u3 / a.u1;
  p = obtain_p(a);

  stress_s = Limit(mid); //径向拉回

  rhos = rho * exp(-3.0 * (stress_s.u1 - s1.u1) /
                   (4.0 * Miu)); //获取弹性极限状态下的密度
  if (abs(rhos - rho) < Eps_rho) {
    ps = p;
    sigma_s =
        -ps + stress_s.u1; //实际密度等于弹性极限密度就认为达到弹性极限状态
    us = u;
    vs = v;
  } else {
    ts = rhos * rho / (rhos - rho);
    ps = (2.0 * Gamma1 * ts * rhos) / (2.0 * ts - Gamma1 * rhos) *
         (const_c * const_c / Gamma1 *
              ((rhos - const_rho) / rhos - (rho - const_rho) / rho) +
          Gamma1 * p / rho);
    sigma_s = -ps + stress_s.u1; //获取弹性极限状态下的物理量
    us = u - sqrt((-p + s1.u1 + ps - stress_s.u1) / ts);
    vs = v;
  }

  U_S.u1 = rhos;
  U_S.u2 = rhos * us; //用原始物理量构造守恒物理量
  U_S.u3 = rhos * vs;
  U_S.u4 = energy(rhos, us, vs, ps);

  return (U_S);
}

/**
 * @brief
 * HLL近似解法器计算x方向右侧弹性极限状态下的物理量（算法跟U_Limit_x_L一样）
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向右侧弹性极限状态下的守恒物理量
 */
inline CVariable U_Limit_x_R(CVariable &a, CVariable &b, CVariable &s1,
                             CVariable &s2) {
  double rho(0), speed(0), ps(0), rhos(0), vs(0), us(0);
  double u(0), v(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double p(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), ts(0);
  double rhor(0), ur(0), pr(0);

  CVariable U_S, stress_s, mid;

  rhor = a.u1;
  ur = a.u2 / a.u1;
  pr = obtain_p(a);

  rho = b.u1;
  u = b.u2 / b.u1;
  v = b.u3 / b.u1;

  p = obtain_p(b);

  sr = max(ur + sonic_speed(a, s1.u1), u + sonic_speed(b, s2.u1));
  sl = -sr;

  U_S = 1.0 / (sr - sl) * (sr * b - sl * a + FU(a, s1) - FU(b, s2));

  mid.u1 = (4.0 / 3.0) * Miu * log(b.u1 / U_S.u1) + s1.u1;
  mid.u2 = -(2.0 / 3.0) * Miu * log(b.u1 / U_S.u1) + s1.u2;
  mid.u3 = s1.u3;
  mid.u4 = 0.0;

  stress_s = Limit(mid);

  rhos = rho * exp(-3.0 * (stress_s.u1 - s1.u1) / (4.0 * Miu));
  if (abs(rhos - rho) < Eps_rho) {
    ps = p;
    sigma_s = -ps + stress_s.u1;
    us = u;
    vs = v;
  } else {
    ts = rhos * rho / (rhos - rho);
    ps = (2.0 * Gamma1 * ts * rhos) / (2.0 * ts - Gamma1 * rhos) *
         (const_c * const_c / Gamma1 *
              ((rhos - const_rho) / rhos - (rho - const_rho) / rho) +
          Gamma1 * p / rho);
    sigma_s = -ps + stress_s.u1;
    us = u + sqrt((-p + s2.u1 + ps - stress_s.u1) / ts);
    vs = v;
  }

  U_S.u1 = rhos;
  U_S.u2 = rhos * us;
  U_S.u3 = rhos * vs;
  U_S.u4 = 0.5 * rhos * us * us + 0.5 * rhos * vs * vs + ps / Gamma1 -
           const_c * const_c * (rhos - const_rho) / Gamma1;

  return (U_S);
}

/**
 * @brief
 * HLLC近似解法器计算x方向右弹性极限状态下的物理量（算法跟HLLC_Limit_x_L一样）
 *
 * @param b i处守恒变量
 * @param mid 预估的右侧星区域的应力
 * @param s2 i处应力
 * CVariable 返回x方向右弹性极限状态下守恒物理量
 */
inline CVariable HLLC_Limit_x_R(CVariable &b, CVariable &mid, CVariable &s2) {
  double rho(0), speed(0), ps(0), rhos(0), vs(0), us(0);
  double u(0), v(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double p(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), ts(0);
  double rhor(0), ur(0), pr(0);

  CVariable U_S, stress_s;

  rho = b.u1;
  u = b.u2 / b.u1;
  v = b.u3 / b.u1;
  p = obtain_p(b);

  stress_s = Limit(mid);

  rhos = rho * exp(-3.0 * (stress_s.u1 - s2.u1) / (4.0 * Miu));
  if (abs(rhos - rho) < Eps_rho) {
    ps = p;
    sigma_s = -ps + stress_s.u1;
    us = u;
    vs = v;
  } else {
    ts = rhos * rho / (rhos - rho);
    ps = (2.0 * Gamma1 * ts * rhos) / (2.0 * ts - Gamma1 * rhos) *
         (const_c * const_c / Gamma1 *
              ((rhos - const_rho) / rhos - (rho - const_rho) / rho) +
          Gamma1 * p / rho);
    sigma_s = -ps + stress_s.u1;
    us = u + sqrt((-p + s2.u1 + ps - stress_s.u1) / ts);
    vs = v;
  }

  U_S.u1 = rhos;
  U_S.u2 = rhos * us;
  U_S.u3 = rhos * vs;
  U_S.u4 = 0.5 * rhos * us * us + 0.5 * rhos * vs * vs + ps / Gamma1 -
           const_c * const_c * (rhos - const_rho) / Gamma1;

  return (U_S);
}

/**
 * @brief HLL近似解法器计算y方向星区域偏应力（算法跟x方向stress_x思路一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向星区域应力
 */
inline CVariable stress_y(CVariable &a, CVariable &b, CVariable &s1,
                          CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0);

  CVariable U_L, U_R, U_S, stress_s, mid;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  sl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  U_S = 1.0 / (sr - sl) * (sr * b - sl * a + GU(a, s1) - GU(b, s2));

  stress_s.u2 = min((4.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u2,
                    (4.0 / 3.0) * Miu * log(b.u1 / U_S.u1) + s2.u2);
  stress_s.u1 = min(-(2.0 / 3.0) * Miu * log(a.u1 / U_S.u1) + s1.u1,
                    -(2.0 / 3.0) * Miu * log(b.u1 / U_S.u1) + s2.u1);

  stress_s.u3 = min(s1.u3, s2.u3);
  stress_s.u4 = 0.0;

  mid = Limit(stress_s);

  return (mid);
}

/**
 * @brief 迎风格式计算x方向数值通量
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向数值通量
 */
inline CVariable X_Upwind(CVariable &a, CVariable &b, CVariable &s1,
                          CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1; //获取原始物理量
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速

  return (0.5 * (FU(a, s1) + FU(b, s1) - 0.5 * sr * (b - a)));
}

/**
 * @brief 迎风格式计算y方向数值通量（算法与X_Upwind一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向数值通量
 */
inline CVariable Y_Upwind(CVariable &a, CVariable &b, CVariable &s1,
                          CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));

  return (0.5 * GU(a, s1) + GU(b, s1) - 0.5 * sr * (b - a));
}

/**
 * @brief L-F格式计算x方向数值通量
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向数值通量
 */
inline CVariable X_LF(CVariable &a, CVariable &b, CVariable &s1,
                      CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1; //获取原始物理量
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //获取波速

  return (0.5 * (FU(a, s1) + FU(b, s1) - sr * (b - a)));
}

/**
 * @brief L-F格式计算y方向数值通量（与X_LF算法一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向数值通量
 */
inline CVariable Y_LF(CVariable &a, CVariable &b, CVariable &s1,
                      CVariable &s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));

  return (0.5 * GU(a, s1) + GU(b, s1) - sr * (b - a));
}

/**
 * @brief HLL近似解法器计算x方向数值通量
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向数值通量
 */
inline CVariable X_HLL(CVariable &a, CVariable &b, CVariable &s1,
                       CVariable &s2) {
  CVariable U_S, stress_L, stress_R;
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1; //获取原始物理量
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速
  sl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  mid1 = stress_x(a, b, s1, s2); //获取星区域应力
  m = vos_limit(mid1); //判断星区域应力是否米塞斯屈服极限
  if ((m - 1.0) > Eps) {

    stress_L = Limit(mid1);
    stress_R = Limit(mid1); //获取弹性极限下的应力

    U_L = U_Limit_x_L(a, b, s1, s2); //获取弹性极限下的守恒量
    U_R = U_Limit_x_R(a, b, s1, s2);

    pl = obtain_p(U_L); //获取弹性极限状态下的压强
    pr = obtain_p(U_R);

    cr = max(U_L.u2 / U_L.u1 + sonic_speed(U_L, stress_L.u1),
             U_R.u2 / U_R.u1 + sonic_speed(U_R, stress_R.u1)); //重新估计波速
    cl = -cr;

    Fflux = 1.0 / (cr - cl) *
            (cr * FU(U_L, stress_L) - cl * FU(U_R, stress_R) +
             cl * cr * (U_R - U_L)); //获取弹性极限状态下的数值通量

    return (Fflux);
  } else {
    Fflux = 1.0 / (sr - sl) *
            (sr * FU(a, mid1) - sl * FU(b, mid1) +
             sl * sr * (b - a)); //获取正常状态下的数值通量
    return (Fflux);
  }
}

/**
 * @brief HLLC近似解法器计算x方向数值通量
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * CVariable 返回x方向数值通量
 */
inline CVariable X_HLLC(CVariable &a, CVariable &b, CVariable &s1,
                        CVariable &s2) {

  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1; //获取原始物理量
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u1;
  sigmar = -pr + s2.u1; //获取总应力

  sl = min(vl - sqrt(Miu / rhol), vr - sqrt(Miu / rhor));
  sr = max(vl + sqrt(Miu / rhol), vr + sqrt(Miu / rhor)); //估计弹性波速

  cr = max(ul + sonic_speed(a, s1.u1),
           ur + sonic_speed(b, s2.u1)); //估计塑性波速
  cl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  Hl = rhol * (cl - ul);
  Hr = rhor * (cr - ur);

  speed_u =
      (-Hl * ul + Hr * ur - sigmal + sigmar) / (Hr - Hl); //获取星区域间断速度
  rhols = Hl / (cl - speed_u);
  rhors = Hr / (cr - speed_u);

  Gl = rhols * (sl - vl);
  Gr = rhors * (sr - vr);
  speed_v = (-Gl * vl + Gr * vr - s1.u3 + s2.u3) / (Gr - Gl);

  sigma_s = sigmar - Hr * (speed_u - ur);  //获取星区域总应力
  stress_xy = s2.u3 - Gr * (speed_v - vr); //获取星区域剪切力

  Els = a.u4 + (speed_u - ul) * (speed_u - sigmal / Hl); //获取星区域能量
  Ers = b.u4 + (speed_u - ur) * (speed_u - sigmar / Hr);

  F_L.u1 = rhols * speed_u;
  F_L.u2 = rhols * speed_u * speed_u - sigma_s;
  F_L.u3 = rhols * speed_u * speed_v - stress_xy; //获取星区域左侧数值通量
  F_L.u4 = (Els - sigma_s) * speed_u - speed_v * stress_xy;

  F_R.u1 = rhors * speed_u;
  F_R.u2 = rhors * speed_u * speed_u - sigma_s;
  F_R.u3 = rhors * speed_u * speed_v - stress_xy; //获取星区域右侧数值通量
  F_R.u4 = (Ers - sigma_s) * speed_u - speed_v * stress_xy;

  if (speed_u >= Eps) { //根据波速选择数值通量
    return (F_L);

  } else {
    return (F_R);
  }
}

/**
 * @brief HLLC近似解法器计算y方向数值通量（算法与X_HLLC一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向数值通量
 */
inline CVariable Y_HLLC(CVariable &a, CVariable &b, CVariable &s1,
                        CVariable &s2) {

  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u2;
  sigmar = -pr + s2.u2;

  sl = min(ul - sqrt(Miu / rhol), ur - sqrt(Miu / rhor));
  sr = max(ul + sqrt(Miu / rhol), ur + sqrt(Miu / rhor));

  cr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  cl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  Hl = rhol * (cl - vl);
  Hr = rhor * (cr - vr);

  speed_v = (-Hl * vl + Hr * vr - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_v);
  rhors = Hr / (cr - speed_v);

  Gl = rhols * (sl - ul);
  Gr = rhors * (sr - ur);
  speed_u = (-Gl * ul + Gr * ur - s1.u3 + s2.u3) / (Gr - Gl);

  sigma_s = sigmar - Hr * (speed_v - vr);
  stress_xy = s2.u3 - Gr * (speed_u - ur);

  Els = a.u4 + (speed_v - vl) * (speed_v - sigmal / Hl);
  Ers = b.u4 + (speed_v - vr) * (speed_v - sigmar / Hr);

  F_L.u1 = rhols * speed_v;
  F_L.u2 = rhols * speed_v * speed_u - stress_xy;
  F_L.u3 = rhols * speed_v * speed_v - sigma_s;
  F_L.u4 = (Els - sigma_s) * speed_v - speed_u * stress_xy;

  F_R.u1 = rhors * speed_v;
  F_R.u2 = rhors * speed_v * speed_u - stress_xy;
  F_R.u3 = rhors * speed_v * speed_v - sigma_s;
  F_R.u4 = (Ers - sigma_s) * speed_v - speed_u * stress_xy;

  if (speed_u >= Eps) {
    return (F_L);

  } else {
    return (F_R);
  }
}

/**
 * @brief 估计X方向波速
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * double 返回X方向波速
 */
inline double X_Uspeed(CVariable &a, CVariable &b, const CVariable s1,
                       const CVariable s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0);

  CVariable U_L, U_R, U_S;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1; //获取原始物理量
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速
  sl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  return (sr);
}

/**
 * @brief HLL近似解法器计算x方向星区域的速度u
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * double 返回x方向星区域的x方向速度u
 */
inline double HLL_X_Ustar(CVariable &a, CVariable &b, CVariable s1,
                          CVariable s2) {
  CVariable U_S, stress_L, stress_R;
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1; //计算原始物理量
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1)); //估计波速
  sl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  stress = stress_x(a, b, s1, s2); //获得星区域应力
  m = vos_limit(stress);           //判断是否超过米塞斯屈服极限
  if ((m - 1.0) > Eps) {

    stress_L = Limit(stress);
    stress_R = Limit(stress); //获取弹性极限下的应力
    U_L = U_Limit_x_L(a, b, s1, s2);
    U_R = U_Limit_x_R(a, b, s1, s2); //获取弹性极限下的守恒物理量

    pl = obtain_p(U_L); //获取弹性极限下的压强
    pr = obtain_p(U_R);

    cr = max(U_L.u2 / U_L.u1 + sonic_speed(U_L, stress_L.u1),
             U_R.u2 / U_R.u1 +
                 sonic_speed(U_R, stress_R.u1)); //估计弹性极限下的波速
    cl = -cr;

    U_S = 1.0 / (cr - cl) *
          (cr * U_R - cl * U_L + FU(U_L, stress_L) -
           FU(U_R, stress_R)); //获取星区域守恒物理量

    return (U_S.u2 / U_S.u1);

  } else {
    U_S = 1.0 / (sr - sl) * (sr * b - sl * a + FU(a, s1) - FU(b, s2));

    return (U_S.u2 / U_S.u1);
  }
}

/**
 * @brief HLLC近似解法器计算x方向星区域的速度u
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * double 返回x方向星区域的速度u
 */
inline double HLLC_X_Ustar(CVariable &a, CVariable &b, const CVariable s1,
                           const CVariable s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1; //计算原始物理量
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u1; //获取总应力
  sigmar = -pr + s2.u1;

  sl = min(vl - sqrt(Miu / rhol), vr - sqrt(Miu / rhor)); //估计弹性波速
  sr = max(vl + sqrt(Miu / rhol), vr + sqrt(Miu / rhor));

  cr = max(ul + sonic_speed(a, s1.u1),
           ur + sonic_speed(b, s2.u1)); //估计塑性波速
  cl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  Hl = rhol * (cl - ul);
  Hr = rhor * (cr - ur);

  speed_u =
      (-Hl * ul + Hr * ur - sigmal + sigmar) / (Hr - Hl); //获取星区域间断速度u
  rhols = Hl / (cl - speed_u);
  rhors = Hr / (cr - speed_u); //获取星区域密度

  Gl = rhols * (sl - vl);
  Gr = rhors * (sr - vr);
  speed_v = (-Gl * vl + Gr * vr - s1.u3 + s2.u3) / (Gr - Gl); //获取星区域密度v

  return (speed_u);
}

/**
 * @brief HLL近似解法器计算y方向数值通量
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * CVariable 返回y方向数值通量
 */
inline CVariable Y_HLL(CVariable &a, CVariable &b, CVariable &s1,
                       CVariable &s2) {
  CVariable U_S, stress_s, mid1, mid2;
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0);

  CVariable U_L, U_R, Gflux;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  sl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  mid1 = stress_y(a, b, s1, s2);

  Gflux =
      1.0 / (sr - sl) * (sr * GU(a, s1) - sl * GU(b, s2) + sl * sr * (b - a));
  return (Gflux);
}

/**
 * @brief 估计y方向波速（算法与X_Uspeed一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * double 返回y方向波速
 */
inline double Y_Vspeed(CVariable &a, CVariable &b, const CVariable s1,
                       const CVariable s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0);

  CVariable U_L, U_R, U_S;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  sl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  return (sr);
}

/**
 * @brief HLL近似解法器计算y方向星区域速度v（算法与HLL_X_Ustar一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * double 返回y方向星区域速度v
 */
inline double HLL_Y_Ustar(CVariable &a, CVariable &b, CVariable s1,
                          CVariable s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0);

  CVariable U_L, U_R, U_S, mid1, mid2;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  sl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  U_S = 1.0 / (sr - sl) * (sr * b - sl * a + GU(a, s1) - GU(b, s2));

  return (U_S.u3 / U_S.u1);
}

/**
 * @brief HLLC近似解法器计算y方向星区域速度v（算法与HLLC_X_Ustar一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * double 返回y方向星区域速度v
 */
inline double HLLC_Y_Ustar(CVariable &a, CVariable &b, const CVariable s1,
                           const CVariable s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u2;
  sigmar = -pr + s2.u2;

  sl = min(ul - sqrt(Miu / rhol), ur - sqrt(Miu / rhor));
  sr = max(ul + sqrt(Miu / rhol), ur + sqrt(Miu / rhor));

  cr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  cl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  Hl = rhol * (cl - vl);
  Hr = rhor * (cr - vr);

  speed_v = (-Hl * vl + Hr * vr - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_v);
  rhors = Hr / (cr - speed_v);

  Gl = rhols * (sl - ul);
  Gr = rhors * (sr - ur);
  speed_u = (-Gl * ul + Gr * ur - s1.u3 + s2.u3) / (Gr - Gl);

  return (speed_v);
}

/**
 * @brief HLL近似解法器计算y方向星区域速度u（算法与HLL_Y_Ustar一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * double 返回y方向星区域速度u
 */
inline double HLL_u_Y_Ustar(CVariable &a, CVariable &b, CVariable s1,
                            CVariable s2) {
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0);

  CVariable U_L, U_R, U_S, mid1, mid2;

  rhol = a.u1;
  rhor = b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  sl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  U_S = 1.0 / (sr - sl) * (sr * b - sl * a + GU(a, s1) - GU(b, s2));

  return (U_S.u2 / U_S.u1);
}

/**
 * @brief HLL近似解法器计算x方向星区域的v（算法与HLL_X_Ustar一致）
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * double 返回x方向星区域的v
 */
inline double HLL_v_X_Ustar(CVariable &a, CVariable &b, CVariable s1,
                            CVariable s2) {
  CVariable U_S, stress_L, stress_R;
  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), m(0), n(0), k(0);

  CVariable U_L, U_R, stress, mid1, mid2;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1));
  sl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  stress = stress_x(a, b, s1, s2);
  m = vos_limit(stress);
  if ((m - 1.0) > Eps) {

    stress_L = Limit(stress);
    stress_R = Limit(stress);

    U_L = U_Limit_x_L(a, b, s1, s2);
    U_R = U_Limit_x_R(a, b, s1, s2);

    pl = obtain_p(U_L);
    pr = obtain_p(U_R);

    cr = max(U_L.u2 / U_L.u1 + sonic_speed(U_L, stress_L.u1),
             U_R.u2 / U_R.u1 + sonic_speed(U_R, stress_R.u1));
    cl = -cr;

    U_S = 1.0 / (cr - cl) *
          (cr * U_R - cl * U_L + FU(U_L, stress_L) - FU(U_R, stress_R));

    return (U_S.u3 / U_S.u1);

  } else {
    U_S = 1.0 / (sr - sl) * (sr * b - sl * a + FU(a, s1) - FU(b, s2));

    return (U_S.u3 / U_S.u1);
  }
}

/**
 * @brief HLLC近似解法器计算y方向星区域速度u（算法与HLLC_Y_Ustar一致）
 *
 * @param a j-1处守恒变量
 * @param b j处守恒变量
 * @param s1 j-1处应力
 * @param s2 j处应力
 * double 返回y方向星区域速度u
 */
inline double HLLC_u_Y_Ustar(CVariable &a, CVariable &b, const CVariable s1,
                             const CVariable s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u2;
  sigmar = -pr + s2.u2;

  sl = min(ul - sqrt(Miu / rhol), ur - sqrt(Miu / rhor));
  sr = max(ul + sqrt(Miu / rhol), ur + sqrt(Miu / rhor));

  cr = max(vl + sonic_speed(a, s1.u2), vr + sonic_speed(b, s2.u2));
  cl = min(vl - sonic_speed(a, s1.u2), vr - sonic_speed(b, s2.u2));

  Hl = rhol * (cl - vl);
  Hr = rhor * (cr - vr);

  speed_v = (-Hl * vl + Hr * vr - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_v);
  rhors = Hr / (cr - speed_v);

  Gl = rhols * (sl - ul);
  Gr = rhors * (sr - ur);
  speed_u = (-Gl * ul + Gr * ur - s1.u3 + s2.u3) / (Gr - Gl);

  return (speed_u);
}

/**
 * @brief HLLC近似解法器计算x方向星区域的v（算法与HLLC_X_Ustar一致）
 *
 * @param a i-1处守恒变量
 * @param b i处守恒变量
 * @param s1 i-1处应力
 * @param s2 i处应力
 * double 返回x方向星区域的v
 */
inline double HLLC_v_X_Ustar(CVariable &a, CVariable &b, const CVariable s1,
                             const CVariable s2) {
  double rhol(0), rhor(0), speed_u(0), speed_v(0), pls(0), prs(0), rhols(0),
      rhors(0);
  double ul(0), ur(0), sl(0), stress_l(0), stress_r(0), cl(0), cr(0), vl(0),
      vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double Hl(0), Hr(0), stress_xy(0), Gl(0), Gr(0);
  double Els(0), Ers(0);

  CVariable F_L, F_R, stress_L, stress_R;

  rhol = a.u1;
  rhor = b.u1;
  ul = a.u2 / a.u1;
  ur = b.u2 / b.u1;
  vl = a.u3 / a.u1;
  vr = b.u3 / b.u1;
  pl = obtain_p(a);
  pr = obtain_p(b);

  sigmal = -pl + s1.u1;
  sigmar = -pr + s2.u1;

  sl = min(vl - sqrt(Miu / rhol), vr - sqrt(Miu / rhor));
  sr = max(vl + sqrt(Miu / rhol), vr + sqrt(Miu / rhor));

  cr = max(ul + sonic_speed(a, s1.u1), ur + sonic_speed(b, s2.u1));
  cl = min(ul - sonic_speed(a, s1.u1), ur - sonic_speed(b, s2.u1));

  Hl = rhol * (cl - ul);
  Hr = rhor * (cr - ur);

  speed_u = (-Hl * ul + Hr * ur - sigmal + sigmar) / (Hr - Hl);
  rhols = Hl / (cl - speed_u);
  rhors = Hr / (cr - speed_u);

  Gl = rhols * (sl - vl);
  Gr = rhors * (sr - vr);
  speed_v = (-Gl * vl + Gr * vr - s1.u3 + s2.u3) / (Gr - Gl);

  return (speed_v);
}

#endif
