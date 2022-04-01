#include<iostream>
#include<string>
using namespace std;

#include"cplex.h"



complex& operator+(complex a, complex b) { return a += b; }
complex operator*(complex a, complex b) { return a *= b; }
int main(){
	complex a{2,1};
	cout << a.imag() << a.real() << endl;
	complex b{1,1};
	complex c(a * b);
	cout << c.imag() <<"\n"<< c.real()<<"\n"<<a.imag()<<a.real();
};