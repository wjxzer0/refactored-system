#include"mana.h"
int main() {
	Employee e{ "brom",1234 };
	Manager mm{ "smith",23456,2 };
	print_list({&e,&mm});
	return 0;
}