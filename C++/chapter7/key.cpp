#include<iostream>
using namespace std;
#include"add.h"
#include"way.cpp"

void print_addr(Address* p) {
	cout << p->name << '\n' << p->number << '\n' << p->street
		<< '\n' << p->town << '\n' << p->state[1] << ' ' << p->zip << endl;
}

void print_addr2(const Address& p) {
	cout << p.name << '\n' << p.number << '\n' << p.street
		<< '\n' << p.town << '\n' << p.state[1] << ' ' << p.zip << endl;
}

void f() {
	Address jd = { "J dad",
61,"South St","New Pro",{'N','J'},"07974" };
	Address JD = jd;
	//print_addr(&jd);
	print_addr2(JD);
	cout << sizeof(JD) << endl;
}

void g() {
	U2 x1;
	U2 x2{ 7 };
	cout << x1.p << endl;
	cout << x2.a << endl;
}

int main() {
	//f();
	g();
	return 0;
}