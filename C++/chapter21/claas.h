#pragma once
#include<iostream>
using namespace std;
class Sate {
public:
	 virtual void center() { cout << "Sate" << endl; }
};

class Disp {
public:
	 virtual void center() { cout << "disp" << endl; }
};

struct DDisp :Disp {
	using Disp::Disp;
	virtual void Disp_cent() = 0;
	void center()override final { Disp_cent(); }
};

struct SSate:Sate
{
	using Sate::Sate;
	virtual void Sate_cent() = 0;
	void center()override final { Sate_cent(); }
};

class comm :public SSate, public DDisp {
public:
	void Disp_cent()override { cout << "disp_center" << endl; }
	void Sate_cent()override { cout << "Sate_center" << endl; }
};