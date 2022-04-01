#pragma once
template<int N>
constexpr int fac() {
	return N * fac<N - 1>();
}

template<>
constexpr int fac<1>() {
	return 1;
}

template<int N>
struct Fac {
	static const int value = N * Fac<N - 1>::value;
};

template<>
struct Fac<1> {
	static const int value = 1;
};
