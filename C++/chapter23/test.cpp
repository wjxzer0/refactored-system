#include<iostream>
#include<fstream>
using namespace std;
int main() {
	ofstream ofs;
	char ch{ '2' };
	ofs.open("um.txt", ios::out);
	ofs << ch;
	ofs.close();
	return 0;
}