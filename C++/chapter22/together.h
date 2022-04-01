#pragma once
#include<iostream>
using namespace std;
class Circle;
class Triangle;

class Shape {
public:
	virtual bool intersect(const Shape&)const = 0;
	virtual bool intersect(const Circle&)const = 0;
	virtual bool intersect(const Triangle&)const = 0;
};

class Circle :public Shape {
public:
	bool intersect(const Shape&)const override;
	virtual bool intersect(const Circle&)const override;
	virtual bool intersect(const Triangle&)const override;
};

class Triangle :public Shape {
public:
	bool intersect(const Shape&)const override;
	virtual bool intersect(const Circle&)const override;
	virtual bool intersect(const Triangle&)const override;
};

bool Circle::intersect(const Shape& s)const { return s.intersect(*this); }
bool Circle::intersect(const Circle&)const { cout << "intersect(Circle,Circle)\n"; return true; }
bool Circle::intersect(const Triangle&)const { cout << "intersect(Circle,Triangle)\n"; return true; }

bool Triangle::intersect(const Shape& s)const { return s.intersect(*this); }
bool Triangle::intersect(const Circle&)const { cout << "intersect(Triangle,circle)\n"; return true; }
bool Triangle::intersect(const Triangle&)const { cout << "intersect(Triangle,Triangle)\n"; return true; }
