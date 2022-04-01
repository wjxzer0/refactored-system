#include"together.h"
#include<vector>
void test(Triangle& t,Circle& c) {
	vector<pair<Shape*, Shape*>>vs{ {&t,&t},{&t,&c},{&c,&c},{&c,&t} };
	for (auto p : vs)
	{
		p.first->intersect(*p.second);
	}
}

int main() {
	Triangle t;
	Circle c;
	cout << typeid(c).name() << endl;
	test(t,c);
	return 0;
}