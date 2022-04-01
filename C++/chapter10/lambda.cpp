#include<iostream>
#include <string>
#include<map>
#include<cctype>
#include<sstream>
using namespace std;
#include "stuctt.h"
#include"way.cpp"




int main(int argc,char* argv[]) {
	switch (argc)
	{
	case 1:
		break;
	case 2:
		ts.set_input(new istringstream{ argv[1] });
		break;
	default:
		error("too many argumengts");
		return 1;
	}
	table["pi"] = 3.1415926535897932385;
	table["e"] = 2.7182818284590452354;

	calculate();

	return no_of_errors;
}