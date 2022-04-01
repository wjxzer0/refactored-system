#include<iostream>
using namespace std;

namespace Chrono {
	enum class Month {jan=1,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec};
	class Date {
	public://�����ӿ�
		class Bad_date {};//�쳣��
		explicit Date(int dd = {}, Month mm = {}, int yy = {});//ѡ��Ĭ��ֵ
		//���޸����ԣ����ڲ���
		int day() const ;
		Month month() const;
		int year()const;

		string string_rep()const;//�ַ�����ʾ
		void char_rep(char s[], int max)const;//C����ַ�����ʾ
		//�޸�����
		Date& add_year(int n);
		Date& add_month(int n);
		Date& add_day(int n);
	private:
		bool is_valid();//���Date�Ƿ��ʾһ������
		int d, m, y;
	};
	bool is_date(int d, Month m, int y);//�Ϸ����ڷ���true
	bool is_leapyear(int y);//���귵��true

	bool operator==(const Date& a, const Date & b);
	bool operator!=(const Date& a, const Date& b);

	const Date& default_date();//Ĭ������

	ostream& operator<<(ostream& os, const Date& d);//��d��ӡ��os
	istream& operator>>(istream& is, Date& d);//��is��ȡDate����d

}
