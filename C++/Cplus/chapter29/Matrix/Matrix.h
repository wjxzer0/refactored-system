#pragma once
#include<iostream>
#include<vector>
#include<cassert>
#include<array>
using namespace std;

#include"Matrix_impl.h"
#include"Matrix_initializer.h"
#include"Matrix_ref.h"


namespace M {
	template<typename T, size_t N>
	class Matrix {
	public:
		static constexpr size_t order = N;
		using value_type = T;
		using iterator = typename std::vector<T>::iterator;
		using const_iterator = typename std::vector<T>::const_iterator;

		Matrix() = default;
		Matrix(Matrix&&) = default;          //移动构造
		Matrix& operator=(Matrix&&) = default;   //移动赋值
		Matrix(Matrix const&) = default;      //拷贝构造
		Matrix& operator=(Matrix const&) = default; //拷贝赋值
		~Matrix() = default;

		template<typename U>
		Matrix(const Matrix_ref<U, N>&);  //构造from Matrix_ref
		template<typename U>
		Matrix& operator=(const Matrix_ref<U, N>&);  //从Matrix_ref赋值

		template<typename... Exts>   //指明每一维大小
		explicit Matrix(Exts... exts);

		Matrix(Matrix_initializer<T, N>);
		Matrix& operator=(Matrix_initializer<T, N>);

		template<typename U>
		Matrix(std::initializer_list<U>) = delete;    //除元素外不使用{}
		template<typename U>
		Matrix& operator=(std::initializer_list<U>) = delete;

		static constexpr size_t get_order() { return N; };   //维数
		size_t extent(const size_t n)const { return desc.extents[n]; }  //第n维元素
		size_t size()const { return elems.size(); }        //元素总数
		const Matrix_slice<N>& descriptor()const { return desc; } //定义下表操作的切片

		T* data() { return elems.data(); }      //平坦元素访问
		const T* data()const { return elems.data(); }



		//下标和切片访问
		template<typename...Args>
		enable_if<Matrix_impl::Requesting_element<Args...>(), T&>//m(i,j,k)用整数进行下标操作
		operator()(Args...args);

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_element<Args...>(),const T&>
			operator()(Args...args)const;

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T,N>>//m(s1,s2,s3)用切片（区间）进行下标操作
			operator()(const Args...args);

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T, N>>
			operator()(const Args...args)const;

		Matrix_ref<T, N - 1>operator[](size_t i) { return row(i); }   //m[i]访问
		Matrix_ref<const T, N - 1>operator[](size_t i)const { return row(i); }

		Matrix_ref<T, N - 1>row(size_t n);//行访问
		Matrix_ref<const T, N - 1>row(size_t n)const;

		Matrix_ref<T, N - 1>col(size_t n);//列访问
		Matrix_ref<const T, N - 1>col(size_t n)const;


		//算术运算
        template<typename F>
		Matrix& apply(F f);  //对每个元素求F（x）

		template<typename F>
		Matrix& apply(const Matrix& m, F f);   //对特定元素执行f(x,m)

		Matrix& operator=(const T& value); //标量赋值

		Matrix& operator+=(const T& value); //标量加
		Matrix& operator-=(const T& value); 
		Matrix& operator*=(const T& value); 
		Matrix& operator/=(const T& value); 
		Matrix& operator%=(const T& value); 

	
		Matrix& operator+=(const Matrix& m);//矩阵加
	
		Matrix& operator-=(const Matrix& m);//矩阵减
		
	private:
		Matrix_slice<N> desc;   //定义N个元素大小的切片
		std::vector<T> elems;        //元素
	};


	template<typename T,size_t N>
	Matrix<T, N> operator+(const Matrix<T, N>& m, const T& val) {
		Matrix<T, N> res = m;
		res += val;
		return res;
	}

	template<typename T, size_t N>
	Matrix<T, N> operator-(const Matrix<T, N>& m, const T& val) {
		Matrix<T, N> res = m;
		res -= val;
		return res;
	}

	template<typename T, size_t N>
	Matrix<T, N> operator*(const Matrix<T, N>& m, const T& val) {
		Matrix<T, N> res = m;
		res *= val;
		return res;
	}

	template<typename T, size_t N>
	Matrix<T, N> operator/(const Matrix<T, N>& m, const T& val) {
		Matrix<T, N> res = m;
		res /= val;
		return res;
	}

	template<typename T, size_t N>
	Matrix<T, N> operator%(const Matrix<T, N>& m, const T& val) {
		Matrix<T, N> res = m;
		res %= val;
		return res;
	}

	template<typename T, size_t N>
	Matrix<T, N> operator+(const Matrix<T, N>& a, const Matrix<T, N>& b) {
		Matrix<T, N> res = a;
		res += b;
		return res;
	}

	template<typename T>
	Matrix<T, 2> operator*(const Matrix<T, 1>& u, const Matrix<T, 1>& v) {
		const size_t n = u.extent(0);
		const size_t m = v.extent(0);
		Matrix<T, 2>res(n, m);
		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j)
				res(i, j) = u[i] * v[j];
		return res;
	}

	template<typename T>
	Matrix<T, 1> operator*(const Matrix<T, 2>& m, const Matrix<T, 1>& v) {
		assert(m.extent(1) == v.extent(0));
		const size_t n = m.extent(0);
		Matrix<T, 1>res(n);
		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != n; ++j)
				res(i) += m[i,j] * v[j];
		return res;
	}

	template<typename T>
	Matrix<T, 2> operator*(const Matrix<T, 2>& m1, const Matrix<T, 2>& m2) {
		const size_t n = m1.extent(0);
		const size_t m = m1.extent(1);
		assert(m == m2.extent(0));
		const size_t p = m2.extent(1);
		Matrix<T, 2>res(n, p);
		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j)
				for (size_t k = 0; k != p; ++k)
					res(i, j) = m1(i, k) * m2(k, j);
		return res;
	}


	template<typename T>
	bool Matrix_type(T& t) {
		//...
		return true;
	}
}


