#include<iostream>
#include <string>
#include<map>
#include<cctype>
#include<sstream>
using namespace std;
enum class Kind :char {
	name,number,end,
	plus='+',minus='-',mul='*',div='/',print=':',assign='=',lp='(',rp=')'
};

struct Token {
	Kind kind;
	string string_value;
	double number_value;
};

class Token_stream {
public:
	Token_stream(istream& s) :ip{ &s }, owns{ false } {cout << "����" << endl; }
	Token_stream(istream* p) :ip{ p }, owns{ true }{cout << "����" << endl; }

	~Token_stream() { close(); cout << "����" << endl; }

	Token get();              //��ȡ��������һ������
	const Token& current();   //�ոն���ĵ���

	void set_input(istream& s) { close(); ip = &s; owns = false; }
	void set_input(istream* p) { close(); ip = p; owns = true; }

private:
	void close() { if (owns) delete ip; }

	istream* ip;
	bool owns;
	Token ct{ Kind::end };
};

