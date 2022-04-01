#include"dd.h"
using namespace Chrono;
//Chrono::Date::Date(int dd, Month mm, int yy) :d{ dd }, m{ mm }, y{yy} {	}
Date& Date::add_month(int n) {
	if (n == 0)return *this;
	if (n > 0) {
		int delay_y = n / 12;
	    m = static_cast<int>(m) + n % 12;
		if (12 < m) {
			++delay_y;
			m -= 12;
		}
		y += delay_y;
		return *this;
	}
	if (n < 0) { cout << "只加不减" << endl; return *this; }
}

Date& Date::add_day(int n) {
	if (n == 0)return *this;
	if (n < 0) { cout << "只加不减" << endl; return *this; }
	if (n > 0) {
		d += n;
		return *this;
	}
}

Date& Date::add_year(int n) {
	if (n < 0) { cout << "只加不减" << endl; return *this; }
	if (n > 0) { y += n; return *this; }
}

int Date::day()const {
	return d;
}

int Date::year()const {
	return y;
}

Month Date::month()const {
	return static_cast<Month>(m);
}



