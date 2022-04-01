#pragma once
#include<iostream>
using namespace std;
class top {
	virtual void write() {};
	virtual void draw() {};

};
class stor :public  top {
	
public:
	int i{};
};

class Trans :public virtual stor {
public:
	void write() { cout << "Trans_write" << endl; }
	void draw() { cout << "Transe_draw" << endl; }
	virtual void run() { cout << "Trans_run" << endl; };
	int j{};
};

class Rece :public virtual stor {
public:
	void write() { cout << "Rece_write" << endl; }
	void draw() {}
	int k{};
};

class Radio :public Trans, public Rece {
public:
	void write()override { Trans::draw(); }
	void draw()override {Rece::write();}
	void run() { cout << "Radio_run" << endl; }
	void laugh() { cout << "haha" << endl; }
	static int m;
};

int Radio::m = 1;