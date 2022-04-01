class complex {
	double re, im;
public:
	complex(double r,double i):re{r},im{i}{}
	complex(double r):re{r},im{0}{}
	complex():re{0},im{0}{}

	double real()const { return re; }
	double imag()const { return im; }
	void real(double d) { re = d; }
	void imag(double d) { im = d; }

	complex& operator+=(complex z) { re += z.re, im += z.im; return *this; }
	complex& operator-=(complex z) { re -= z.re, im -= z.im; return *this; }
	//complex& operator*=(complex);
	complex& operator/=(complex);


};

complex operator*=(complex a, complex b) {
	double re = a.real() * b.real() - a.imag() * b.imag();
	double im = a.real() * b.imag() + a.imag() * b.real();
	complex r(re, im);
	return r;
};
