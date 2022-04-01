#pragma once
#include"cont.h"
#include<iostream>
#include<list>
#include<typeinfo>
using namespace std;

namespace vct {

	class Container {
	public:
		virtual double& operator[](int) = 0;
		virtual int size()const = 0;
		virtual ~Container() {}

	};

	class Vector {
	private:
		double* elem;
		int sz;
	public:
		Vector(int s) :elem{ new double[s] }, sz{ s }{
			for (int i = 0; i != s; ++i)
				elem[i] = 0;
			cout << "����v" << endl;
		}
		Vector(initializer_list<double>);
		double& operator[](int i) { return elem[i]; }
		~Vector() {
			cout << "ɾ��" << elem[1] << endl;
			delete[]elem;
		}
		int size()const { return sz; }
	};

	class Vector_container :public Container {
		Vector v;
	public:
		Vector_container(int s) :v(s) {
			cout << "����V" << endl;
		}
		~Vector_container() {
			cout << "����V";
		}

		double& operator[](int i) { return v[i]; }
		int size()const { return v.size(); }
	};


	class List_container :public Container {
		list<double>id;
	public:
		List_container() {}
		List_container(initializer_list<double>il) :id{ il } {
			cout << "����L" << endl;
		}
		~List_container() {}
		double& operator[](int i);
		int size()const { return id.size(); }
	};
	
	void use(Container& c);
	void g();

	void h();

	void f();
	
}