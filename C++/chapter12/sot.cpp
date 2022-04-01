#include<iostream>
#include<vector>
using namespace std;
#include"Us.h"
#include<algorithm>
vector<User>heads = {
	{"Risi","dmr",11271},
	{"Seth","ravi",11272},
	{"Szyman","tgs",11273},
	{"Schr","nls",11274},
	{"Schr","nls",11275},
	{"Kerni","bwk",11276}
};

void print_id(vector<User>& v) {
	for (auto& x : v)
		cout << x.name << '\t' << x.id << '\t' << x.dept << endl;
}

int cmp1(const void* p, const void* q) {
	return strcmp(static_cast<const User*>(p)->name, static_cast<const User*>(q)->name);
}

int cmp2(const void* p, const void* q) {
	return static_cast<const User*>(p)->dept, static_cast<const User*>(q)->dept;
}

using CFT = int(const void*, const void*);

void ssort(void* base, size_t n, size_t sz, CFT cmp) {
	for (int gap = n / 2; 0 < gap; gap /= 2)
		for (int i = gap; i != n; i++)
			for (int j = i - gap; 0 <= j; j -= gap) {
				char* b = static_cast<char*>(base);
				char* pj = b + j * sz;
				char* pjg = b + (j + gap) * sz;
				if (cmp(pjg, pj) < 0) {
					for (int k = 0; k != sz; k++) {
						char temp = pj[k];
						pj[k] = pjg[k];
						pjg[k] = temp;
					}
				}
			}				
}

void ss() {
	cout << "Head in alpha oder\n";
	//ssort(heads, 6, sizeof(User), cmp1);
	print_id(heads);
	cout << "Head in dept oder\n";
	//ssort(heads, 6, sizeof(User), cmp2);
	print_id(heads);
}

void s() {
	cout << "Head in alpha oder\n";
	sort(heads.begin(), heads.end(), [](const User& x, const User& y) {return x.name < y.name; });
	print_id(heads);
	cout << "Head in dept oder\n";
	sort(heads.begin(), heads.end(), [](const User& x, const User& y) {return x.dept < y.dept; });
	print_id(heads);
}

int main() {
	print_id(heads);
	//ss();
	s();
	return 0;
}