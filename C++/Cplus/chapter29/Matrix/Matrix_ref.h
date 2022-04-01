#pragma once
using namespace std;

namespace M {
template<typename T,size_t N>
class Matrix_ref {
public:
	Matrix_ref(const Matrix_slice<N>& s, T* p) :desc{ s }, ptr{ p }{}//ºÜÏñMatrix
private:
	Matrix_slice<N>  desc;
	T* ptr;
};
}