#include<iostream>
#include"Matrix.h"
using namespace std;


void test01() {
	Matrix<double, 2>m1{ {1,2},{3,4} };
	cout << m1.get_order() << endl;
	//cout << m1.extent(1) << endl;
	//cout << m1.size() << endl;
	//Matrix<double, 1>m2 = m1[1];
	//cout << m2.size() << endl;
}

int main() {}