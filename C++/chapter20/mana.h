#pragma once
#include<iostream>
#include<list>
using namespace std;
class Employee {
public:
	Employee(const string& name, int dept) :family_name{ name }, department{ dept }{}
	virtual void print() { cout << family_name << '\t' << department << endl; }

private:
	string first_name, family_name;
	int department;
};

class Manager :public Employee {
public:
	Manager(const string& name, int dept, int lvl) :level{ lvl }, Employee(name, dept){ }
	void print() { Employee::print(); cout << "\tlevel" << level << endl; }
private:
	list<Employee>group;
	int level;
};

void print_list(const list<Employee*>& s) {
	for (auto x : s)
		x->print();
}