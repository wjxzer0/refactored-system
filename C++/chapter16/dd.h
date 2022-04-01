#include<iostream>
using namespace std;

namespace Chrono {
	enum class Month {jan=1,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec};
	class Date {
	public://公共接口
		class Bad_date {};//异常类
		explicit Date(int dd = {}, Month mm = {}, int yy = {});//选择默认值
		//非修改属性，用于查找
		int day() const ;
		Month month() const;
		int year()const;

		string string_rep()const;//字符串表示
		void char_rep(char s[], int max)const;//C风格字符串表示
		//修改属性
		Date& add_year(int n);
		Date& add_month(int n);
		Date& add_day(int n);
	private:
		bool is_valid();//检查Date是否表示一个日期
		int d, m, y;
	};
	bool is_date(int d, Month m, int y);//合法日期返回true
	bool is_leapyear(int y);//闰年返回true

	bool operator==(const Date& a, const Date & b);
	bool operator!=(const Date& a, const Date& b);

	const Date& default_date();//默认日期

	ostream& operator<<(ostream& os, const Date& d);//将d打印到os
	istream& operator>>(istream& is, Date& d);//从is读取Date存入d

}
