#pragma once

#include "Stdafx.hpp"

//using namespace std;

class Cell {

public:
  Cell(int ii, int jj) : i(ii), j(jj){}

//   ~Cell()=default;

private:
  int i;
  int j;

public:

  int getcell_i() const { return (i); };

  int getcell_j() const { return (j); };

};
