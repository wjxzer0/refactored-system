#include "vc.h"

double& vct::List_container::operator[](int i) {
	for (auto& x : id) {
		if (i == 0)return x;
		--i;
	}
	throw out_of_range("List container");
}

vct::Vector::Vector(std::initializer_list<double> lst) :elem{ new double[lst.size()] }, sz{ static_cast<int>(lst.size()) }{
	copy(lst.begin(), lst.end(), elem);
	cout << "¹¹Ôì" << typeid(5).name() << std::endl;
}

void vct::use(Container& c)
{
	const int sz = c.size();
	for (int i = 0; i != sz; ++i)
		cout << c[i] << "\n";
}
void vct::g() {
	Vector_container vc(9);
	use(vc);
}

void vct::h() {
	List_container lc = { 1,2,3,4,5,6,7,8,9 };
	cout << lc[2] << "test" << endl;
	use(lc);
}

void vct::f() {

}