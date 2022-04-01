#include "Variable.h"

CVariable::CVariable() : u1(0), u2(0), u3(0), u4(0) {}

CVariable::CVariable(double a, double b, double c, double d) {
  u1 = a;
  u2 = b;
  u3 = c;
  u4 = d;
}

CVariable &CVariable::operator=(const CVariable &v) {
  u1 = v.u1;
  u2 = v.u2;
  u3 = v.u3;
  u4 = v.u4;

  return (*this);
}

CVariable CVariable::operator+(const CVariable &v) {
  return CVariable(u1 + v.u1, u2 + v.u2, u3 + v.u3, u4 + v.u4);
}

CVariable CVariable::operator-(const CVariable &v) {
  return CVariable(u1 - v.u1, u2 - v.u2, u3 - v.u3, u4 - v.u4);
}

CVariable operator*(double c, const CVariable &v) {
  return CVariable(c * v.u1, c * v.u2, c * v.u3, c * v.u4);
}

double operator*(const CVariable &v1, const CVariable &v2) {
  return (v1.u1 * v2.u1 + v1.u2 * v2.u2 + v1.u3 * v2.u3 + v1.u4 * v2.u4);
}

double CVariable::GetAbsMax() {

  double s = fabs(u1);

  if (fabs(u2) > s)
    s = fabs(u2);
  if (fabs(u3) > s)
    s = fabs(u3);
  if (fabs(u4) > s)
    s = fabs(u4);

  return s;
}
