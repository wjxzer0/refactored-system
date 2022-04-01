#include"strr.h"
char* str_ing::expand(const char* ptr, int n) {
	char* p = new char[n];
	strcpy(p,ptr);
	return p;
}

void str_ing::String::copy_from(const String& x) {
	if (x.sz <= short_max) {
		memcpy(this, &x, sizeof(x));
		ptr = ch;
	}
	else {
		ptr = expand(x.ptr, x.sz + 1);
		sz = x.sz;
		space = 0;
	}
}

void str_ing::String::move_from(String& x) {
	if (x.sz <= short_max) {
		memcpy(this, &x, sizeof(x));
		ptr = ch;
	}
	else {
		ptr = x.ptr;
		sz = x.sz;
		space = x.space;
		x.ptr = x.ch;    //x=""
		x.sz = 0;
		x.ch[0] = 0;
	}
}