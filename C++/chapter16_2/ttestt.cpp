#include<iostream>
using namespace std;

class Date {
	int d;
public:
	Date(int dd) { d = dd ? dd : 10; }
	void prt() { cout << d << endl; }
};



int main() {
	Date da(1);
	da.prt();
	return 0;
}