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
	Token_stream(istream& s) :ip{ &s }, owns{ false } {cout << "构造" << endl; }
	Token_stream(istream* p) :ip{ p }, owns{ true }{cout << "构造" << endl; }

	~Token_stream() { close(); cout << "析构" << endl; }

	Token get();              //读取并返回下一个单次
	const Token& current();   //刚刚读入的单次

	void set_input(istream& s) { close(); ip = &s; owns = false; }
	void set_input(istream* p) { close(); ip = p; owns = true; }

private:
	void close() { if (owns) delete ip; }

	istream* ip;
	bool owns;
	Token ct{ Kind::end };
};

