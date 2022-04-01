#include"claas.h"
#include"over.h"
/*
void g(Disp* p) {
	cout << sizeof(p) << '\t' << sizeof(*p) << endl;
}

void f(Sate* q) {
	cout << sizeof(q) << '\t' << sizeof(*q) << endl;
}
*/
void g(Trans* t);
void test01() {
	cout << sizeof(Disp) << '\t' << sizeof(Sate) << '\t' << sizeof(comm) << endl;
	comm pd;
	pd.DDisp::center();
	pd.SSate::center();
	pd.Disp_cent();
	pd.Sate_cent();
	//g(new comm);
	//f(&pd);
}

void test02() {
	cout <<sizeof(stor)<<'\t'<< sizeof(Trans)<<'\t'<< sizeof(Rece) << '\t' << sizeof(Radio) << endl;
	Radio r;
	r.draw();
	r.write();
	g(new Radio);
	g(new Trans);
}

void g(Trans* t) { t->run(); t->write(); }


void test03() {
	Radio r;
	top* pw = &r;
	cout << typeid(pw).name() << endl;
	cout << typeid(*pw).name() << endl;
	stor* pd = dynamic_cast<stor*>(pw);
	cout << typeid(pd).name() << endl;
	Trans t;
	//Radio& pw = dynamic_cast<Radio&>(t);
	if (pd == nullptr) { cout << "¿Õ" << endl; return; }
	cout << "·Ç¿Õ" << endl;
}

int main() {
	cout << "¿ªÊ¼" << endl;
	//test01();
	//test02();
	test03();
	return 0;
}