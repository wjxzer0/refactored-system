#include "Stdafx.hpp"
#include "TwoDimElasticPlastic.h"

using namespace std;

int main() {
  //初值向量，分别包含rho , u ,v, p
  vector<double> Ls, Rs;

  //初值向量，分别包含s_xx,s_yy,s_xy
  vector<double> Lstress, Rstress;

  Ls.resize(4);
  Rs.resize(4);

  Lstress.resize(3);
  Rstress.resize(3);

  Ls[0] = 2.7;
  Ls[1] = 2.0;
  Ls[2] = 0.0;
  Ls[3] = 1.0;

  Rs[0] = 2.7;
  Rs[1] = 0.0;
  Rs[2] = 0.0;
  Rs[3] = 1.0;

  Lstress[0] = 0.0;
  Lstress[1] = 0.0;
  Lstress[2] = 0.0;

  Rstress[0] = 0.0;
  Rstress[1] = 0.0;
  Rstress[2] = 0.0;

  CTwoDimElasticPlastic test(-2.0, 2.0, -2.0, 2.0, 400, 400, 0.005, Ls, Rs,
                             Lstress, Rstress, 3000);

  test.MainSolve();

  test.Output();

  return 0;
}