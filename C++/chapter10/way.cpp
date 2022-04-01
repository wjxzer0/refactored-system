#include"stuctt.h"
int no_of_errors{ 0 };
Token_stream ts{ cin };
map <string, double>table;
double prim(bool get);
double term(bool get);

Token_stream ct;

double error(const string& s) {
	no_of_errors++;
	cerr << "error:" << s << '\n';
	return 1;
}

double expr(bool get)//�Ӽ� 
{
	double left = term(get);

	for (;;) {//forever
		switch (ts.current().kind)
		{
		case Kind::plus:
			left += term(true);
			break;
		case Kind::minus:
			left -= term(true);
		default:
			return left;

		}
	}
}

double term(bool get)//�˳� 
{
	double left = term(get);

	for (;;) {//forever
		switch (ts.current().kind)
		{
		case Kind::mul:
			left *= prim(true);
			break;
		case Kind::div:
			if (auto d = prim(true)) {
				left /= d;
				break;
			}
			return error("divide by 0");
		default:
			return left;

		}
	}
}

double prim(bool get) {//���������

	if (get) ts.get();//��ȡ��һ������

	switch (ts.current().kind)
	{
	case Kind::number://����������
	{double v = ts.current().number_value;
	ts.get();
	return v; }
	case Kind::name:
	{double& v = table[ts.current().string_value];//�ҵ���Ӧ��
	if (ts.get().kind == Kind::assign)v = expr(true);//������=�������
	return v;
	}
	case Kind::minus:				//һԪ����
		return -prim(true);
	case Kind::lp:
	{auto e = expr(true);
	if (ts.current().kind != Kind::rp)return error("')'expected");
	ts.get();					//�Ե��ˡ�����
	return e;
	}
	default:
		return error("primary expected");
	}
}

void calculate() {
	for (;;) {
		ts.get();
		if (ts.current().kind == Kind::end)break;
		if (ts.current().kind == Kind::print)continue;
		cout << expr(false) << '\n';
	}
}


Token Token_stream::get() {
	char ch;
	*ip >> ch;
	switch (ch)
	{
	case 0:
		return ct = { Kind::end };
	case ';':
	case '*':
	case '/':
	case '+':
	case '-':
	case '(':
	case ')':
	case '=':
		return ct = { static_cast<Kind>(ch) };
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	case '.':
		ip->putback(ch);
		*ip >> ct.number_value;
		ct.kind = Kind::number;
		return ct;
	default:
		if (isalpha(ch)) {
			ip->putback(ch);
			*ip >> ct.string_value;
			ct.kind = Kind::name;
			return ct;
		}
	}
}

