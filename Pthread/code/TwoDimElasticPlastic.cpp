#include "TwoDimElasticPlastic.h"
#include "Stdafx.hpp"

CTwoDimElasticPlastic::CTwoDimElasticPlastic(){};

CTwoDimElasticPlastic::~CTwoDimElasticPlastic(){};

CTwoDimElasticPlastic::CTwoDimElasticPlastic(
    double left, double right, double bottom, double top, int xNum, int yNum,
    double time, vector<double> &leftstate, vector<double> &rightstate,
    vector<double> &leftstress, vector<double> &rightstress, double y_e) {

  m_left = left;
  m_right = right;
  m_top = top;
  m_bottom = bottom;

  m_xNum = xNum;
  m_yNum = yNum;

  m_xH = (right - left) / xNum;
  m_yH = (top - bottom) / yNum;

  m_Time = time;

  xEnd = m_xNum + xBegin;
  yEnd = m_yNum + yBegin;

  xTotal = xBegin + xEnd;
  yTotal = yBegin + yEnd;

  m_ini_lrho = leftstate[0];
  m_ini_lu = leftstate[1];
  m_ini_lv = leftstate[2];
  m_ini_lp = leftstate[3];

  m_ini_rrho = rightstate[0];
  m_ini_ru = rightstate[1];
  m_ini_rv = rightstate[2];
  m_ini_rp = rightstate[3];

  m_ini_lsx = leftstress[0];
  m_ini_lsy = leftstress[1];
  m_ini_lsxy = leftstress[2];

  m_ini_rsx = rightstress[0];
  m_ini_rsy = rightstress[1];
  m_ini_rsxy = rightstress[2];

  m_ini_Y = y_e;

  Grid.Resize(xTotal, yTotal);

  U0.Resize(xTotal, yTotal);
  S0.Resize(xTotal, yTotal);
  Y_E.Resize(xTotal, yTotal);
  tempY.Resize(xTotal, yTotal);

  tempU0.Resize(xTotal, yTotal);
  tempU1.Resize(xTotal, yTotal);
  tempU2.Resize(xTotal, yTotal);
  tempU3.Resize(xTotal, yTotal);

  tempS0.Resize(xTotal, yTotal);
  tempS1.Resize(xTotal, yTotal);
  tempS2.Resize(xTotal, yTotal);
  tempS3.Resize(xTotal, yTotal);

  tempY0.Resize(xTotal, yTotal);
  tempY1.Resize(xTotal, yTotal);
  tempY2.Resize(xTotal, yTotal);
  tempY3.Resize(xTotal, yTotal);

  FU0.Resize(xTotal, yTotal);
  GU0.Resize(xTotal, yTotal);

  Fflux.Resize(xTotal + 1, yTotal);
  X_speed.Resize(xTotal + 1, yTotal);
  v_X_speed.Resize(xTotal + 1, yTotal);

  Gflux.Resize(xTotal, yTotal + 1);
  Y_speed.Resize(xTotal, yTotal + 1);
  u_Y_speed.Resize(xTotal, yTotal + 1);

  distance_phi.Resize(xTotal, yTotal);
  distance_phi_temp.Resize(xTotal, yTotal);

  // left_U0.Resize(xTotal, yTotal);
  // right_U0.Resize(xTotal, yTotal);
  // top_U0.Resize(xTotal, yTotal);
  // bottom_U0.Resize(xTotal, yTotal);

  // left_S0.Resize(xTotal, yTotal);
  // right_S0.Resize(xTotal, yTotal);
  // top_S0.Resize(xTotal, yTotal);
  // bottom_S0.Resize(xTotal, yTotal);

  X_dericition_delta_phi.Resize(xTotal, yTotal);
  Y_dericition_delta_phi.Resize(xTotal, yTotal);

  // Grid information

  Grid(xBegin, yBegin).x = m_left + 0.5 * m_xH;
  Grid(xBegin, yBegin).y = m_bottom + 0.5 * m_yH;

  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      Grid(i, j).x = Grid(xBegin, yBegin).x + (i - xBegin) * m_xH;
      Grid(i, j).y = Grid(xBegin, yBegin).y + (j - yBegin) * m_yH;
    }
  }
}

void CTwoDimElasticPlastic::Initial1() {

  //初始间断
  double x0 = 0.0;
  double coff = 0.0;

  //初始距离函数
  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      distance_phi(i, j) = distance_function(Grid(i, j).x, Grid(i, j).y);
    }
  }

  //初始赋值
  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      if (distance_phi(i, j) < Eps) {
        if (Grid(i, j).x < (x0 + coff * Grid(i, j).y)) {
          U0(i, j).u1 = m_ini_lrho;
          U0(i, j).u2 = m_ini_lrho * m_ini_lu;
          U0(i, j).u3 = m_ini_lrho * m_ini_lv;
          U0(i, j).u4 = energy(m_ini_lrho, m_ini_lu, m_ini_lv,
                               m_ini_lp); //刚性气体状态方程

          S0(i, j).u1 = m_ini_lsx;
          S0(i, j).u2 = m_ini_lsy;
          S0(i, j).u3 = m_ini_lsxy;
          S0(i, j).u4 = 0.0;

          Y_E(i, j).u1 = m_ini_Y;
          Y_E(i, j).u2 = 0.0;
          Y_E(i, j).u3 = 0.0;
          Y_E(i, j).u4 = 0.0;
        }

        else {
          U0(i, j).u1 = m_ini_rrho;
          U0(i, j).u2 = m_ini_rrho * m_ini_ru;
          U0(i, j).u3 = m_ini_rrho * m_ini_rv;
          U0(i, j).u4 = energy(m_ini_rrho, m_ini_ru, m_ini_rv,
                               m_ini_rp); //刚性气体状态方程

          S0(i, j).u1 = m_ini_rsx;
          S0(i, j).u2 = m_ini_rsy;
          S0(i, j).u3 = m_ini_rsxy;
          S0(i, j).u4 = 0.0;

          Y_E(i, j).u1 = m_ini_Y;
          Y_E(i, j).u2 = 0.0;
          Y_E(i, j).u3 = 0.0;
          Y_E(i, j).u4 = 0.0;
        }
      } else {
        U0(i, j).u1 = 2.7;
        U0(i, j).u2 = 0.0;
        U0(i, j).u3 = 0.0;
        U0(i, j).u4 = energy(2.7, 0.0, 0.0, 1.0); //刚性气体状态方程

        S0(i, j).u1 = 0.0;
        S0(i, j).u2 = 0.0;
        S0(i, j).u3 = 0.0;
        S0(i, j).u4 = 0.0;

        Y_E(i, j).u1 = m_ini_Y;
        Y_E(i, j).u2 = 0.0;
        Y_E(i, j).u3 = 0.0;
        Y_E(i, j).u4 = 0.0;
      }
    }
  }
}

void CTwoDimElasticPlastic::Initial2() {

  //初始赋值
  for (int i = xBegin + 200; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      {
        U0(i, j).u1 = m_ini_rrho;
        U0(i, j).u2 = m_ini_rrho * m_ini_ru;
        U0(i, j).u3 = m_ini_rrho * m_ini_rv;
        U0(i, j).u4 =
            energy(m_ini_rrho, m_ini_ru, m_ini_rv, m_ini_rp); //刚性气体状态方程

        S0(i, j).u1 = m_ini_rsx;
        S0(i, j).u2 = m_ini_rsy;
        S0(i, j).u3 = m_ini_rsxy;
        S0(i, j).u4 = 0.0;

        Y_E(i, j).u1 = m_ini_Y;
        Y_E(i, j).u2 = 0.0;
        Y_E(i, j).u3 = 0.0;
        Y_E(i, j).u4 = 0.0;
      }
    }
  }
  for (int i = xBegin; i != xBegin + 200; i++) {
    for (int j = yBegin + 100; j != yEnd - 100; j++) {
      {
        U0(i, j).u1 = m_ini_lrho;
        U0(i, j).u2 = m_ini_lrho * m_ini_lu;
        U0(i, j).u3 = m_ini_lrho * m_ini_lv;
        U0(i, j).u4 =
            energy(m_ini_lrho, m_ini_lu, m_ini_lv, m_ini_lp); //刚性气体状态方程

        S0(i, j).u1 = m_ini_lsx;
        S0(i, j).u2 = m_ini_lsy;
        S0(i, j).u3 = m_ini_lsxy;
        S0(i, j).u4 = 0.0;

        Y_E(i, j).u1 = m_ini_Y;
        Y_E(i, j).u2 = 0.0;
        Y_E(i, j).u3 = 0.0;
        Y_E(i, j).u4 = 0.0;
      }
    }
  }
  for (int i = xBegin; i != xBegin + 200; i++) {
    for (int j = yBegin; j != yBegin + 100; j++) {
      {
        U0(i, j).u1 = m_ini_rrho;
        U0(i, j).u2 = m_ini_rrho * m_ini_ru;
        U0(i, j).u3 = m_ini_rrho * m_ini_rv;
        U0(i, j).u4 =
            energy(m_ini_rrho, m_ini_ru, m_ini_rv, m_ini_rp); //刚性气体状态方程

        S0(i, j).u1 = m_ini_rsx;
        S0(i, j).u2 = m_ini_rsy;
        S0(i, j).u3 = m_ini_rsxy;
        S0(i, j).u4 = 0.0;

        Y_E(i, j).u1 = m_ini_Y;
        Y_E(i, j).u2 = 0.0;
        Y_E(i, j).u3 = 0.0;
        Y_E(i, j).u4 = 0.0;
      }
    }
  }
  for (int i = xBegin; i != xBegin + 200; i++) {
    for (int j = yEnd - 100; j != yEnd; j++) {
      {
        U0(i, j).u1 = m_ini_rrho;
        U0(i, j).u2 = m_ini_rrho * m_ini_ru;
        U0(i, j).u3 = m_ini_rrho * m_ini_rv;
        U0(i, j).u4 =
            energy(m_ini_rrho, m_ini_ru, m_ini_rv, m_ini_rp); //刚性气体状态方程

        S0(i, j).u1 = m_ini_rsx;
        S0(i, j).u2 = m_ini_rsy;
        S0(i, j).u3 = m_ini_rsxy;
        S0(i, j).u4 = 0.0;

        Y_E(i, j).u1 = m_ini_Y;
        Y_E(i, j).u2 = 0.0;
        Y_E(i, j).u3 = 0.0;
        Y_E(i, j).u4 = 0.0;
      }
    }
  }
}

void CTwoDimElasticPlastic::SetBoundary(double time) {
  // left-> inflow condition
  double rho(0), u(0), v(0), p(0), a(0), E(0), H(0);

  for (int i = 0; i != xBegin; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      U0(i, j) = U0(xBegin, j);

      S0(i, j) = S0(xBegin, j);

      Y_E(i, j) = S0(xBegin, j);
    }
  }
  // right-> outflow condition
  for (int j = yBegin; j != yEnd; j++) {
    U0(xEnd, j) = U0(xEnd - 1, j);
    U0(xEnd + 1, j) = U0(xEnd - 2, j);
    U0(xEnd + 2, j) = U0(xEnd - 3, j);

    S0(xEnd, j) = S0(xEnd - 1, j);
    S0(xEnd + 1, j) = S0(xEnd - 2, j);
    S0(xEnd + 2, j) = S0(xEnd - 3, j);

    Y_E(xEnd, j) = S0(xEnd - 1, j);
    Y_E(xEnd + 1, j) = S0(xEnd - 2, j);
    Y_E(xEnd + 2, j) = S0(xEnd - 3, j);
  }
  // bottom
  for (int i = xBegin; i != xEnd; i++) {
    for (int j = 0; j != yBegin; j++)

      U0(i, 2) = U0(i, yBegin);
    U0(i, 1) = U0(i, yBegin + 1);
    U0(i, 0) = U0(i, yBegin + 2);

    S0(i, 2) = S0(i, yBegin);
    S0(i, 1) = S0(i, yBegin + 1);
    S0(i, 0) = S0(i, yBegin + 2);

    Y_E(i, 2) = S0(i, yBegin);
    Y_E(i, 1) = S0(i, yBegin + 1);
    Y_E(i, 0) = S0(i, yBegin + 2);
  }

  // top
  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yEnd; j != yTotal; j++) {
      U0(i, j) = U0(i, 2.0 * yEnd - j - 1);

      S0(i, j) = S0(i, 2.0 * yEnd - j - 1);

      Y_E(i, j) = S0(i, 2.0 * yEnd - j - 1);
    }
  }
}

void CTwoDimElasticPlastic::GetX_Uspeed() {
  for (int j = yBegin; j != (yEnd); j++) {
    for (int i_flux = xBegin; i_flux != (xEnd + 1); i_flux++) {
      X_speed(i_flux, j).u1 = HLLC_X_Ustar(U0(i_flux - 1, j), U0(i_flux, j),
                                           S0(i_flux - 1, j), S0(i_flux, j));
      X_speed(i_flux, j).u2 = 0.0;
      X_speed(i_flux, j).u3 = 0.0;
      X_speed(i_flux, j).u4 = 0.0;

      v_X_speed(i_flux, j).u1 = HLLC_v_X_Ustar(
          U0(i_flux - 1, j), U0(i_flux, j), S0(i_flux - 1, j), S0(i_flux, j));
      v_X_speed(i_flux, j).u2 = 0.0;
      v_X_speed(i_flux, j).u3 = 0.0;
      v_X_speed(i_flux, j).u4 = 0.0;
    }
  }
}

void CTwoDimElasticPlastic::GetX_flux() {
  for (int j = yBegin; j != yEnd; j++) {
    for (int i_flux = xBegin; i_flux != (xEnd + 1); i_flux++) {
      Fflux(i_flux, j) = X_HLLC(U0(i_flux - 1, j), U0(i_flux, j),
                                S0(i_flux - 1, j), S0(i_flux, j));
    }
  }
}

void CTwoDimElasticPlastic::GetY_flux() {
  for (int i = xBegin; i != xEnd; i++) {
    for (int j_flux = yBegin; j_flux != (yEnd + 1); j_flux++) {
      Gflux(i, j_flux) = Y_HLLC(U0(i, j_flux - 1), U0(i, j_flux),
                                S0(i, j_flux - 1), S0(i, j_flux));
    }
  }
}

void CTwoDimElasticPlastic::GetY_Uspeed() {

  for (int i = xBegin; i != (xEnd); i++) {
    for (int j_flux = yBegin; j_flux != (yEnd + 1); j_flux++) {

      Y_speed(i, j_flux).u1 = HLLC_Y_Ustar(U0(i, j_flux - 1), U0(i, j_flux),
                                           S0(i, j_flux - 1), S0(i, j_flux));
      Y_speed(i, j_flux).u2 = 0.0;
      Y_speed(i, j_flux).u3 = 0.0;
      Y_speed(i, j_flux).u4 = 0.0;

      u_Y_speed(i, j_flux).u1 = HLLC_u_Y_Ustar(
          U0(i, j_flux - 1), U0(i, j_flux), S0(i, j_flux - 1), S0(i, j_flux));
      u_Y_speed(i, j_flux).u2 = 0.0;
      u_Y_speed(i, j_flux).u3 = 0.0;
      u_Y_speed(i, j_flux).u4 = 0.0;
    }
  }
}

// void CTwoDimElasticPlastic::Real_area_index(double time, int x_point, int
// point_y) {

//   for (int i = x_point-5; i != x_point+5; i++) {
//     for (int j = point_y-5; j != point_y+5; j++) {

//       if (distance_phi(i, j) < 0.0) {
//         Index.push_back(Cell(i, j));
//       }
//     }
//   }
// }

void CTwoDimElasticPlastic::Virtual_assignment(double time) {
  vector<Cell> Index;
  Center interface;
  CVariable U_temp, S_temp;
  Clevel_set_method test(-2.0, 2.0, -2.0, 2.0, 400, 400, time, 0.00005, U0,
                         distance_phi);

  //  CMatrix<double> a,b;
  test.Initial();
  // test.X_dericition_delta_phi()=a;


cout<<12345<<endl;

 test.X_dericition_delta_phi();
  test.Y_dericition_delta_phi();

  cout<<1<<endl;
  //  X_dericition_delta_phi = test.X_dericition_delta_phi();
  //  Y_dericition_delta_phi = test.Y_dericition_delta_phi();

   cout<<X_dericition_delta_phi(5,5)<<endl;

  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {

      if (0.0 < distance_phi(i, j) < 1.5 * (m_xH + m_yH)) {
        Index.push_back(Cell(i, j));
      }
    }
  }
  for (auto &b : Index) {
    int m = b.getcell_i();
    int n = b.getcell_j();

    int p = 0, q = 0;

    interface = interface_points(
        Grid(m, n).x, Grid(m, n).y, X_dericition_delta_phi(m, n),
        Y_dericition_delta_phi(m, n), distance_phi(m, n));

    int k = floor((interface.x - Grid(xBegin, yBegin).x) / m_xH + xBegin);
    int l = floor((interface.y - Grid(xBegin, yBegin).y) / m_yH + yBegin);

    double max =
        abs((Grid(k, l).x - Grid(m, n).x) * X_dericition_delta_phi(m, n) +
            (Grid(k, l).y - Grid(m, n).y) * Y_dericition_delta_phi(m, n));

    for (int i = k - 3; i != k + 3; i++) {
      for (int j = l - 3; j != l + 3; j++) {
        if (distance_phi(i, j) < 0.0) {

          if (max <
              abs((Grid(i, j).x - Grid(m, n).x) * X_dericition_delta_phi(m, n) +
                  (Grid(i, j).y - Grid(m, n).y) *
                      Y_dericition_delta_phi(m, n))) {
            max = abs(
                (Grid(i, j).x - Grid(m, n).x) * X_dericition_delta_phi(m, n) +
                (Grid(i, j).y - Grid(m, n).y) * Y_dericition_delta_phi(m, n));
            p = i;
            q = j;
          }
        }
      }
    }
    U_temp =
        single_Riemann_solver(U0(p, q), S0(p, q), X_dericition_delta_phi(p, q),
                              Y_dericition_delta_phi(p, q));
    U0(m, n) = U_temp;

    S0(m, n).u1 = 0.0;
    S0(m, n).u2 = 0.0;
    S0(m, n).u3 = 0.0;
    S0(m, n).u4 = 0.0;
  }

  test.TempSolver(3);
  //  distance_phi_temp = test.level_set_Output();
  //  distance_phi = distance_phi_temp;
}

double CTwoDimElasticPlastic::GetTimeStep(void) {
  double rho, u, v, p, a;
  double maxtemp(0), temp(0);

  for (int i = xBegin; i != xEnd; i++) {
    for (int j = yBegin; j != yEnd; j++) {
      rho = U0(i, j).u1;
      u = U0(i, j).u2 / U0(i, j).u1;
      v = U0(i, j).u3 / U0(i, j).u1;
      p = obtain_p(U0(i, j));

      a = sonic_speed(U0(i, j), S0(i, j).u1);

      temp = (fabs(u) + a) / m_xH + (fabs(v) + a) / m_yH;

      if (temp > maxtemp) {
        maxtemp = temp;
      }
    }
  }

  return (CFL / maxtemp);
}

void CTwoDimElasticPlastic::Timestep(double t) {
  //边界赋值
  SetBoundary(t);
  double m, n;

  //通量赋值
  for (int i = 0; i != xTotal; i++) {
    for (int j = 0; j != yTotal; j++) {
      FU0(i, j) = FU(U0(i, j), S0(i, j));
      GU0(i, j) = GU(U0(i, j), S0(i, j));
    }
  }

  GetX_flux();

  GetX_Uspeed();

  GetY_flux();

  GetY_Uspeed();
}

void CTwoDimElasticPlastic::Rungekutta(double t, double m_tstep) {
  int i, j;
  double w1, w2, w3, w4, b, n;

  for (i = xBegin; i != xEnd; i++) {
    for (j = yBegin; j != yEnd; j++) {
      tempU0(i, j) = U0(i, j);
      tempS0(i, j) = S0(i, j);
      tempY0(i, j) = Y_E(i, j); //赋值进行计算
    }
  }

  Timestep(t);

  for (i = xBegin; i != xEnd + 1; i++) {
    for (j = yBegin; j != yEnd + 1; j++) {
      tempU1(i, j) =
          tempU0(i, j) +
          m_tstep * (-1.0 / m_xH * (Fflux(i + 1, j) - Fflux(i, j)) -
                     1.0 / m_yH * (Gflux(i, j + 1) - Gflux(i, j))); //第一阶RK
    }
  }

  ///////////////////////////////////
  for (i = xBegin; i != xEnd; i++) {
    for (j = yBegin; j != yEnd; j++) {
      U0(i, j) = tempU1(i, j);
    }
  }

  Timestep(t + m_tstep);

  for (i = xBegin; i != xEnd + 1; i++) {
    for (j = yBegin; j != yEnd + 1; j++) {
      tempU2(i, j) =
          0.75 * tempU0(i, j) + 0.25 * tempU1(i, j) +
          0.25 * m_tstep *
              (-1.0 / m_xH * (Fflux(i + 1, j) - Fflux(i, j)) -
               1.0 / m_yH * (Gflux(i, j + 1) - Gflux(i, j))); //第二阶RK
    }
  }

  ///////////////////////////////////////
  for (i = xBegin; i != xEnd; i++) {
    for (j = yBegin; j != yEnd; j++) {
      U0(i, j) = tempU2(i, j);
    }
  }

  Timestep(t + 0.5 * m_tstep);

  for (i = xBegin; i != xEnd + 1; i++) {
    for (j = yBegin; j != yEnd + 1; j++) {
      tempU3(i, j) =
          1.0 / 3.0 * tempU0(i, j) + 2.0 / 3.0 * tempU2(i, j) +
          2.0 / 3.0 * m_tstep *
              (-1.0 / m_xH * (Fflux(i + 1, j) - Fflux(i, j)) -
               1.0 / m_yH * (Gflux(i, j + 1) - Gflux(i, j))); //第三阶RK
    }
  }

  for (i = xBegin; i != xEnd + 1; i++) {
    for (j = yBegin; j != yEnd + 1; j++) {

      tempY(i, j).u1 =
          sqrt((pow(abs(tempS0(i, j).u1), 2.0)) + //计算线性硬化系数
               (pow(abs(tempS0(i, j).u2), 2.0)));
      tempY(i, j).u2 = 0.0;
      tempY(i, j).u3 = 0.0;
      tempY(i, j).u4 = 0.0;

      if (tempY(i, j).u1 <=
          sqrt(2.0 / 3.0) * tempY0(i, j).u1) { //判断是否发生线性硬化

        tempS1(i, j).u1 =
            tempS0(i, j).u1 +
            2.0 * Miu *
                (2.0 / 3.0 * (m_tstep / m_xH) *
                     (X_speed(i + 1, j).u1 - X_speed(i, j).u1) -
                 1.0 / 3.0 * (m_tstep / m_yH) *
                     (Y_speed(i, j + 1).u1 - Y_speed(i, j).u1)); //更新应力
        tempS1(i, j).u2 = tempS0(i, j).u2 +
                          2.0 * Miu *
                              (2.0 / 3.0 * (m_tstep / m_yH) *
                                   (Y_speed(i, j + 1).u1 - Y_speed(i, j).u1) -
                               1.0 / 3.0 * (m_tstep / m_xH) *
                                   (X_speed(i + 1, j).u1 - X_speed(i, j).u1));
        tempS1(i, j).u3 =
            tempS0(i, j).u3 +
            Miu * m_tstep *
                ((v_X_speed(i + 1, j).u1 - v_X_speed(i, j).u1) / m_xH +
                 (u_Y_speed(i, j + 1).u1 - u_Y_speed(i, j).u1) / m_yH);

        tempS1(i, j).u4 = 0.0;

      } else if (tempY(i, j).u1 <= sqrt(2.0 / 3.0) * Y0) {

        tempS1(i, j).u1 = tempS0(i, j).u1 +
                          2.0 * Miu_p *
                              (2.0 / 3.0 * (m_tstep / m_xH) *
                                   (X_speed(i + 1, j).u1 - X_speed(i, j).u1) -
                               1.0 / 3.0 * (m_tstep / m_yH) *
                                   (Y_speed(i, j + 1).u1 - Y_speed(i, j).u1));
        tempS1(i, j).u2 = tempS0(i, j).u2 +
                          2.0 * Miu_p *
                              (2.0 / 3.0 * (m_tstep / m_yH) *
                                   (Y_speed(i, j + 1).u1 - Y_speed(i, j).u1) -
                               1.0 / 3.0 * (m_tstep / m_xH) *
                                   (X_speed(i + 1, j).u1 - X_speed(i, j).u1));
        tempS1(i, j).u3 =
            tempS0(i, j).u3 +
            Miu_p * m_tstep *
                ((v_X_speed(i + 1, j).u1 - v_X_speed(i, j).u1) / m_xH +
                 (u_Y_speed(i, j + 1).u1 - u_Y_speed(i, j).u1) / m_yH);

        tempS1(i, j).u4 = 0.0;

      } else {
        tempS1(i, j) = Limit(tempS0(i, j));
      }
    }
  }

  for (i = xBegin; i != xEnd; i++) {
    for (j = yBegin; j != yEnd; j++) {
      n = vos_limit(tempS1(i, j));
      if ((n - 1.0) > Eps) {
        S0(i, j) = Limit(tempS1(i, j));
        U0(i, j) = tempU3(i, j);

      } else {
        S0(i, j) = tempS1(i, j);
        U0(i, j) = tempU3(i, j);
      }
    }
  }

  for (i = xBegin; i != xEnd; i++) {
    for (j = yBegin; j != yEnd; j++) {
      Y_E(i, j).u1 =
          max(tempY0(i, j).u1,
              sqrt(3.0 / 2.0) *
                  sqrt((pow(abs(tempS1(i, j).u1), 2.0)) +
                       (pow(abs(tempS1(i, j).u2), 2.0)))); //更新线性硬化系数
      Y_E(i, j).u2 = 0.0;
      Y_E(i, j).u3 = 0.0;
      Y_E(i, j).u4 = 0.0;
    }
  }
}

void CTwoDimElasticPlastic::MainSolve() {
  double m_tstep(0);
  double sum(0);

  Initial1();

  while (sum < m_Time - 1.0e-6) {
    m_tstep = GetTimeStep();
    Virtual_assignment(m_tstep);
    Rungekutta(sum, m_tstep);

    sum = sum + m_tstep;

    cout << ceil(sum / m_tstep) << ", "
         << "Time: " << sum << "s" << endl;
  }
}

void CTwoDimElasticPlastic::Output() {
  int k = 200;
  double rho, u, v, p, e, S_x, S_y, S_xy, sigma_x, sigma_y, Y, m;

  ofstream out_final("result.plt", ios::out);
  out_final.setf(ios::fixed | ios::showpoint);
  out_final << "Variables=\"x\",\"y\",\"rho\",\"u\",\"v\",\"p\",\"S_x\",\"S_"
               "y\",\"S_xy\",\"Sigma_x\",\"Sigma_y\",\"Y\" "
            << endl;
  out_final << "Zone I=" << m_xNum << "  j=" << m_yNum << "  f=point" << endl;

  for (int j = yBegin; j != yEnd; j++) {
    for (int i = xBegin; i != xEnd; i++) {
      m = vos_limit(S0(i, j));
      if ((m - 1.0) > Eps) {
        rho = U0(i, j).u1;
        u = U0(i, j).u2 / U0(i, j).u1;
        v = U0(i, j).u3 / U0(i, j).u1;
        p = obtain_p(U0(i, j));

        S_x = Limit(S0(i, j)).u1;
        S_y = Limit(S0(i, j)).u2;
        S_xy = Limit(S0(i, j)).u3;

        sigma_x = -p + S_x;
        sigma_y = -p + S_y;
        Y = Y_E(i, j).u1;
      } else {
        rho = U0(i, j).u1;
        u = U0(i, j).u2 / U0(i, j).u1;
        v = U0(i, j).u3 / U0(i, j).u1;
        p = obtain_p(U0(i, j));

        S_x = S0(i, j).u1;
        S_y = S0(i, j).u2;
        S_xy = S0(i, j).u3;

        sigma_x = -p + S_x;
        sigma_y = -p + S_y;
        Y = Y_E(i, j).u1;
      }

      out_final << Grid(i, j).x << " " << Grid(i, j).y << " " << rho << " " << u
                << " " << v << " " << p << " " << S_x << " " << S_y << " "
                << S_xy << " " << sigma_x << " " << sigma_y << " " << Y << endl;
    }
  }
  out_final.close();

  ofstream out_final1("density.dat", ios::out);
  out_final1.setf(ios::fixed | ios::showpoint);
  out_final1 << "Variables=\"x\",\"rho\" " << endl;
  out_final1 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {
    rho = U0(i, k).u1;
    out_final1 << Grid(i, k).x << " " << rho << " " << endl;
  }

  out_final1.close();

  ofstream out_final2("velocity_X.dat", ios::out);
  out_final2.setf(ios::fixed | ios::showpoint);
  out_final2 << "Variables=\"x\",\"u\" " << endl;
  out_final2 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    u = U0(i, k).u2 / U0(i, k).u1;

    out_final2 << Grid(i, k).x << " " << u << " " << endl;
  }

  out_final2.close();

  ofstream out_final3("pressure.dat", ios::out);
  out_final3.setf(ios::fixed | ios::showpoint);
  out_final3 << "Variables=\"x\",\"p\" " << endl;
  out_final3 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    p = Gamma1 * U0(i, k).u4 -
        0.5 * Gamma1 *
            (U0(i, k).u2 * U0(i, k).u2 / U0(i, k).u1 +
             U0(i, k).u3 * U0(i, k).u3 / U0(i, k).u1) +
        const_c * const_c * (U0(i, k).u1 - const_rho);

    out_final3 << Grid(i, k).x << " " << p << " " << endl;
  }

  out_final3.close();

  ofstream out_final4("velocity_y.dat", ios::out);
  out_final4.setf(ios::fixed | ios::showpoint);
  out_final4 << "Variables=\"x\",\"v\" " << endl;
  out_final4 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    v = U0(i, k).u3 / U0(i, k).u1;

    out_final4 << Grid(i, k).x << " " << v << " " << endl;
  }

  out_final4.close();

  ofstream out_final5("stress_x.dat", ios::out);
  out_final5.setf(ios::fixed | ios::showpoint);
  out_final5 << "Variables=\"x\",\"S_x\"" << endl;
  out_final5 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    m = vos_limit(S0(i, k));
    if ((m - 1.0) > Eps) {

      S_x = Limit(S0(i, k)).u1;

    } else {

      S_x = S0(i, k).u1;
    }

    out_final5 << Grid(i, k).x << " " << S_x << " " << endl;
  }

  out_final5.close();

  ofstream out_final6("stress_y.dat", ios::out);
  out_final6.setf(ios::fixed | ios::showpoint);
  out_final6 << "Variables=\"x\",\"S_y\"" << endl;
  out_final6 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {
    m = vos_limit(S0(i, k));
    if ((m - 1.0) > Eps) {

      S_y = Limit(S0(i, k)).u2;

    } else {

      S_y = S0(i, k).u2;
    }

    out_final6 << Grid(i, k).x << " " << S_y << " " << endl;
  }

  out_final6.close();

  ofstream out_final7("stress_xy.dat", ios::out);
  out_final7.setf(ios::fixed | ios::showpoint);
  out_final7 << "Variables=\"x\",\"S_xy\"" << endl;
  out_final7 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    m = vos_limit(S0(i, k));
    if ((m - 1.0) > Eps) {

      S_xy = Limit(S0(i, k)).u3;

    } else {

      S_xy = S0(i, k).u3;
    }

    out_final7 << Grid(i, k).x << " " << S_xy << " " << endl;
  }

  out_final7.close();

  ofstream out_final8("sigma_x.dat", ios::out);
  out_final8.setf(ios::fixed | ios::showpoint);
  out_final8 << "Variables=\"x\",\"Sigma_x\" " << endl;
  out_final8 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    m = vos_limit(S0(i, k));
    if ((m - 1.0) > Eps) {

      p = Gamma1 * U0(i, k).u4 -
          0.5 * Gamma1 *
              (U0(i, k).u2 * U0(i, k).u2 / U0(i, k).u1 +
               U0(i, k).u3 * U0(i, k).u3 / U0(i, k).u1) +
          const_c * const_c * (U0(i, k).u1 - const_rho);

      S_x = Limit(S0(i, k)).u1;
      sigma_x = -p + S_x;

    } else {
      p = Gamma1 * U0(i, k).u4 -
          0.5 * Gamma1 *
              (U0(i, k).u2 * U0(i, k).u2 / U0(i, k).u1 +
               U0(i, k).u3 * U0(i, k).u3 / U0(i, k).u1) +
          const_c * const_c * (U0(i, k).u1 - const_rho);

      S_x = S0(i, k).u1;
      sigma_x = -p + S_x;
    }

    out_final8 << Grid(i, k).x << " " << sigma_x << " " << endl;
  }

  out_final8.close();

  ofstream out_final9("sigma_y.dat", ios::out);
  out_final9.setf(ios::fixed | ios::showpoint);
  out_final9 << "Variables=\"x\",\"Sigma_y\" " << endl;
  out_final9 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    m = vos_limit(S0(i, k));
    if ((m - 1.0) > Eps) {
      p = Gamma1 * U0(i, k).u4 -
          0.5 * Gamma1 *
              (U0(i, k).u2 * U0(i, k).u2 / U0(i, k).u1 +
               U0(i, k).u3 * U0(i, k).u3 / U0(i, k).u1) +
          const_c * const_c * (U0(i, k).u1 - const_rho);

      S_y = Limit(S0(i, k)).u2;
      sigma_y = -p + S_y;
    } else {
      p = Gamma1 * U0(i, k).u4 -
          0.5 * Gamma1 *
              (U0(i, k).u2 * U0(i, k).u2 / U0(i, k).u1 +
               U0(i, k).u3 * U0(i, k).u3 / U0(i, k).u1) +
          const_c * const_c * (U0(i, k).u1 - const_rho);

      S_y = S0(i, k).u2;
      sigma_y = -p + S_y;
    }

    out_final9 << Grid(i, k).x << " " << sigma_y << " " << endl;
  }

  out_final9.close();

  ofstream out_final10("limit.dat", ios::out);
  out_final10.setf(ios::fixed | ios::showpoint);
  out_final10 << "Variables=\"x\",\"limit\" " << endl;
  out_final10 << "Zone I=" << m_xNum << "  f=point" << endl;

  for (int i = xBegin; i != xEnd; i++) {

    m = vos_limit(S0(i, k));

    out_final10 << Grid(i, k).x << " " << m << " " << endl;
  }

  out_final10.close();
}
