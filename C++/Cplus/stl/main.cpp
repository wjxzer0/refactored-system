#include<iostream>
#include<vector>
#include<algorithm>
#include<functional>
using namespace std;

template<typename T>
void show(const vector<T> m) {
	for (auto x : m) {
		cout << x << '\t';
	}
}

template<typename T>
void pr(const T val) {
	cout << val<<'\t';
}

template<typename T>
class cmp{
public:
	bool operator()(const T a, const T b)const;
};

template<typename T>
bool cmp<T>::operator()(const T a,const T b)const {
	return a > b ? true : false;
}

void test01() {
	vector<int> m = { 1,2,3,7,6,5 };
	m.push_back(10);
	m.resize(12, 2);
	for_each(m.begin(), m.end(), pr<int>);
	cout << endl;
	sort(m.begin(), m.end(), greater<int>{});
	for (auto x : m) {
		cout << x << '\t';
	}
}

void test02() {
	vector<int> m{1,3,5,2,7,4};
	show(m);
	vector<int>m1(&m[1],&m[4]);
	for_each(m1.begin(), m1.end(),pr<int>);
	cout << endl << typeid(m.begin()).name() << endl << typeid(&m).name() << endl;
	
}

void test03() {
	vector<int> m{ 1,3,5,2,7,4 };
	cout << m.front() << endl;
	cout << m.back() << endl;
	cout << m.at(5);
}

int main() {
	test01();
	//test02();
	//test03();
	return 0;
}

