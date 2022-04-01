#pragma once
#include<iostream>
using namespace std;
void intval() {
	for (char c; cin >> c;)
		cout << "the value of'" << c << "'is " << int{ c } << "\n";
}

void digits() {
	for (int i = 0; i != 10; ++i)
		cout << static_cast<char>('0' + i);
}

void f()
{
	char v[] = "animal";
	char* p = v;
	cout << strlen(p) << endl;
	cout << strlen(v) << endl;
}

void test() {
	/*
		signed char sc = -160;
		unsigned char uc = sc;
		std::cout << uc << std::endl;
		std::cout << "\10";*/


	int a[]{ 1,2,3,4,5 };
	cout << a[4] << endl;
	int* p = &a[0];
	cout << *p << "\n" << a[0] << "\n" << ++ * p << endl;
	cout << *p << "\n" << a << "\n" << ++p << endl;

	/*
	char p[] = "abcd";
	std::cout << p << std::endl;
	p[3] = 'z';
	std::cout << sizeof(p) <<p<< std::endl;*/
	/*
	std::string s = "\\w\\\\w";
	std::string s1 = R"(\w\\w)";
	std::cout << s1 << std::endl;*/
}

int ma[3][5];



void init_ma() {
	for (int i = 0; i != 3; ++i)
		for (int j = 0; j != 5; j++)
			ma[i][j] = 10 * i + j;
}

void print_ma() {
	for (int i = 0; i != 3; ++i) {
		for (int j = 0; j != 5; j++)
			cout << ma[i][j] << "\t";
		cout << "\n";
	}
}
