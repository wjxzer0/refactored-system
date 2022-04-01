/********************************************
// File name: Global.hpp

//这个模块是定义全局变量，包括各种材料常数

// Current version: 1.0
// Author: Zeng zhiqiang(1257371357@qq.com)
// Date: September,21 2020
***********************************************/

#ifndef _GLOBAL_H_
#define _GLOBAL_H_
#pragma once

#include "Matrix.hpp"
#include "Variable.h"
#include <vector>

using namespace std;

//材料比热比
const double Gamma = 2.67;

// Gamma - 1.0;
const double Gamma1 = 1.67; // Gamma - 1.0;

// CFL数
const double CFL = 0.4;

//弹性模量
const double Kpa = 7.4e5;

//弹性杨氏模量
const double Miu = 2.65e5;

//塑性杨氏模量
const double Miu_p = 2.65e5;

//给定的可容许误差
const double Eps = 1.0e-14;

//给定的可容许误差
const double Weno_Eps = 1.0e-6;

//用于x方向边界赋值的网格
const int xBegin = 3;

//用于y方向边界赋值的网格
const int yBegin = 3;

//给定初始声速
const double const_c = 538;

//给定初始密度
const double const_rho = 2.71;

//塑性屈服极限
const double Y0 = 3000;

//给定的可容许密度误差
const double Eps_rho = 1e-3;

//pai值
const double pai=3.141592653;

#endif