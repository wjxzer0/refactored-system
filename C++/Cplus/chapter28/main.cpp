#include<iostream>
#include"contral.h"
#include"plus.h"
using namespace std;

template<bool B, typename T, typename N>
using Conditional = typename conditional<B, T, N>::type;

void test() {
	cout << conditional<true, Squart, Cube>::type{}(99) << endl;
	cout << Conditional<false, Squart, Cube>{}(99) << endl;
	cout << ((2 > 1) ? Squart{}(99) : Cube{}(99)) << endl;
}

void test01() {
	constexpr int x = fac<5>();
	constexpr int x2 = Fac<7>::value;
	cout << x <<'\t'<<x2 << endl;
}

void test02(int* n) {
	n[0] = 1;
}

int main() {
	//test01();
	int* n = new int[10]{};
	cout << n[0]<<"\t";
	test02(n);
	cout << n[0] << endl;
	return 0;
}