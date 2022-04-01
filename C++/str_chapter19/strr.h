#include<iostream>
using namespace std;
namespace str_ing {
	class String {
	public:
		String();
		explicit String(const char* p);
		String(const String&);
		String& operator=(const String&);
		String(String&& x);
		String& operator=(String&& x);
		~String() { if (short_max < sz)delete[] ptr; }

		//操作
		char& operator[](int n) { return ptr[n]; }   //未经检查的访问元素
		char operator[](int n)const { return ptr[n]; }

		char& at(int n) { check(n); return ptr[n]; }   //执行边界检查的访问元素
		char at(int n) const{ check(n); return ptr[n]; }

		String& operator+=(char c);        //在末尾添加c

		const char* c_ptr() { return ptr; }  //C风格字符串访问
		const char* c_ptr() const{ return ptr; }

		int size() const { return sz; }
		int capacity() const { return (sz <= short_max) ? short_max : sz + space; }

	private:
		static const int short_max = 15;
		int sz;
		char* ptr;
		union { 
			int space;
			char ch[short_max+1];//结尾的0
		};
		void check(int n) const //边界检查
		{
			if (n < 0 || sz <= n)
				throw std::out_of_range("String::at()");
		}
		void copy_from(const String& x);
		void move_from(String& x);
	};

	char* expand(const char* ptr, int n);

	String::String() :sz{ 0 }, ptr{ch} {
		ch[0] = 0;
	}

	String::String(const char* p) : sz{ static_cast<int>(strlen(p)) },
		ptr{ (sz <= short_max) ? ch : new char[sz + 1] }, 
		space{0} {
		strcpy(ptr, p);
	}

	String::String(const String& x) {
		copy_from(x);
	}

	String::String(String&& x) {
		move_from(x);
	}

	String& String::operator=(const String& x) {
		if (this == &x)return *this;
		char* p = (short_max < sz) ? ptr : 0;
		copy_from(x);
		delete[] p;
		return *this;
	}

	String& String::operator=(String&& x) {
		if (this == &x)return *this;
		if (short_max < sz)delete[]ptr;
		move_from(x);
		return *this;
	}

	String& String::operator+=(char c) {
		if (sz == short_max) {
			int n = sz + sz + 2;
			ptr = expand(ptr, n);
			space = n - sz - 2;
		}
		else if (short_max < sz) {
			if (space == 0) {
				int n = sz + sz + 2;
				char* p = expand(ptr, n);
				delete[]ptr;
				ptr = p;
				space = n - sz - 2;
			}
			else
				--space;
		}
		ptr[sz] = c;
		ptr[++sz] = 0;
		return *this;
	}

	ostream& operator<<(ostream& os, const String& s);
	istream& operator>>(istream& is, const String& s);
	bool operator==(const String& a, const String& b);
	bool operator!=(const String& a, const String& b);
	char* begin(String& x);
	char* end(String& x);
	const char* begin(const String& x);
	const char* end(const String& x);
	String& operator+=(String& a, const String& b);
	String operator+(const String& a, const String& b);

}