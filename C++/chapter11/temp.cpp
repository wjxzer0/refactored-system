#include<iostream>
#include<string>
using namespace std;

constexpr int isqrt_help(int sq, int d, int a) {
	return sq <= a ? isqrt_help(sq + d, d + 2, a) : d;
}

constexpr int isqrt(int x) {
	return isqrt_help(1, 3, x)/2-1;
}

void f() {
	for (int c; cin >> c;)
		cout << isqrt(c) << endl;
}

constexpr const int* addr(const int &r) { return &r; }

void g() {
	static const int x{5};
	constexpr const int* y = addr(x);
	//constexpr const int yy = *y;
	cout << y << "\n" << &x << endl;
}

template<class T, int N>void h(T(&r)[N]) {};

int main() {

	//f();
	//g();
	using intt = int;
	intt a[]{1,2,3,4};
	h(a);
	const char* day[]{ "one","twe","three" };
	cout << sizeof(day)/sizeof(&day[0]) << endl;
	cout << numeric_limits<int>::max() << endl;
	int* p = new int{ 10 };
	int* q = p;
	cout << *p << *q << endl;
	cout << p << q << endl;
	delete p;
	cout << p << '\t' << *q << endl;

	return 0;
}