#include "Approximate_Riemann_solver.hpp"
#include "Global.hpp"
#include "Variable.h"
#include "cell.hpp"

using namespace std;

inline double distance_function(const double x, const double y) {
  return ((pow(x, 2)) + pow(y, 2) - 1.0);
}

inline double X_Normal_projection(double ul, double vl,
                                  double X_dericition_delta_phi,
                                  double Y_dericition_delta_phi) {

  double m(0), n(0), ur(0), vr(0);

  m = X_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  n = Y_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  ur = ul * m + vl * n;

  return (ur);
}

inline double Y_Normal_projection(double ul, double vl,
                                  double X_dericition_delta_phi,
                                  double Y_dericition_delta_phi) {

  double m(0), n(0), ur(0), vr(0);

  m = X_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  n = Y_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  vr = -ul * n + vl * m;

  return (vr);
}
inline double X_Normal_Projection_back(double ul, double vl,
                                       double X_dericition_delta_phi,
                                       double Y_dericition_delta_phi) {

  double m(0), n(0), ur(0), vr(0);

  m = X_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  n = Y_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  ur = ul * m - vl * n;

  return (ur);
}

inline double Y_Normal_Projection_back(double ul, double vl,
                                       double X_dericition_delta_phi,
                                       double Y_dericition_delta_phi) {

  double m(0), n(0), ur(0), vr(0);

  m = X_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  n = Y_dericition_delta_phi /
      sqrt(pow(X_dericition_delta_phi, 2.0) + pow(Y_dericition_delta_phi, 2.0));

  vr = ul * n + vl * m;

  return (vr);
}

inline CVariable single_Riemann_solver(CVariable &a, CVariable &s1,
                                       double X_dericition_delta_phi,
                                       double Y_dericition_delta_phi) {

  double rhol(0), rhor(0), speed(0), pls(0), prs(0), rhols(0), rhors(0);
  double ul(0), ur(0), sl(0), cl(0), cr(0), vl(0), vr(0);
  double pl(0), pr(0), sr(0), sigmal(0), sigmar(0), sigma_s(0);
  double el(0), er(0), us(0), vs(0);

  CVariable U_L, U_R, stress, mid1, mid2, Fflux;

  rhol = a.u1;
  ul = X_Normal_projection(a.u2 / a.u1, a.u3 / a.u1, X_dericition_delta_phi,
                           Y_dericition_delta_phi);
  vl = Y_Normal_projection(a.u2 / a.u1, a.u3 / a.u1, X_dericition_delta_phi,
                           Y_dericition_delta_phi);
  pl = obtain_p(a);
  sigmal = -pl + s1.u1;
  el = internal_energy(rhol, pl);

  cl = sonic_speed(a, s1.u1);

  rhor = 1.0 / (1.0 / rhol + (-1.0 - sigmal) / (pow(rhol * cl, 2)));
  er = el + 0.5 * ((1.0 / rhol - 1.0 / rhor) * (sigmal + -1.0));

  if (abs(sigmal) > 1.0) {

    ur = ul + sqrt((sigmal + 1.0) * (1.0 / rhol - 1.0 / rhor));

  } else {
    ur = ul - sqrt((sigmal + 1.0) * (1.0 / rhol - 1.0 / rhor));
  }

  us = X_Normal_Projection_back(ur, vr, X_dericition_delta_phi,
                                Y_dericition_delta_phi);

  U_R.u1 = rhor;
  U_R.u2 = rhor * us;
  U_R.u3 = rhor * vl;
  U_R.u4 = er + 0.5 * rhor * us * us + 0.5 * rhor * vl * vl;

  return (U_R);
}

inline Center interface_points(double x, double y,
                               double X_dericition_delta_phi,
                               double Y_dericition_delta_phi,
                               double distance_phi) {

  double x_0, y_0;
  Center Interface;

  Interface.x = x - X_dericition_delta_phi * abs(distance_phi);
  Interface.y = y - Y_dericition_delta_phi * abs(distance_phi);

  return (Interface);
}