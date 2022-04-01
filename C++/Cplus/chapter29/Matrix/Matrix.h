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
		Matrix(Matrix&&) = default;          //�ƶ�����
		Matrix& operator=(Matrix&&) = default;   //�ƶ���ֵ
		Matrix(Matrix const&) = default;      //��������
		Matrix& operator=(Matrix const&) = default; //������ֵ
		~Matrix() = default;

		template<typename U>
		Matrix(const Matrix_ref<U, N>&);  //����from Matrix_ref
		template<typename U>
		Matrix& operator=(const Matrix_ref<U, N>&);  //��Matrix_ref��ֵ

		template<typename... Exts>   //ָ��ÿһά��С
		explicit Matrix(Exts... exts);

		Matrix(Matrix_initializer<T, N>);
		Matrix& operator=(Matrix_initializer<T, N>);

		template<typename U>
		Matrix(std::initializer_list<U>) = delete;    //��Ԫ���ⲻʹ��{}
		template<typename U>
		Matrix& operator=(std::initializer_list<U>) = delete;

		static constexpr size_t get_order() { return N; };   //ά��
		size_t extent(const size_t n)const { return desc.extents[n]; }  //��nάԪ��
		size_t size()const { return elems.size(); }        //Ԫ������
		const Matrix_slice<N>& descriptor()const { return desc; } //�����±��������Ƭ

		T* data() { return elems.data(); }      //ƽ̹Ԫ�ط���
		const T* data()const { return elems.data(); }



		//�±����Ƭ����
		template<typename...Args>
		enable_if<Matrix_impl::Requesting_element<Args...>(), T&>//m(i,j,k)�����������±����
		operator()(Args...args);

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_element<Args...>(),const T&>
			operator()(Args...args)const;

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T,N>>//m(s1,s2,s3)����Ƭ�����䣩�����±����
			operator()(const Args...args);

		template<typename...Args>
		enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T, N>>
			operator()(const Args...args)const;

		Matrix_ref<T, N - 1>operator[](size_t i) { return row(i); }   //m[i]����
		Matrix_ref<const T, N - 1>operator[](size_t i)const { return row(i); }

		Matrix_ref<T, N - 1>row(size_t n);//�з���
		Matrix_ref<const T, N - 1>row(size_t n)const;

		Matrix_ref<T, N - 1>col(size_t n);//�з���
		Matrix_ref<const T, N - 1>col(size_t n)const;


		//��������
        template<typename F>
		Matrix& apply(F f);  //��ÿ��Ԫ����F��x��

		template<typename F>
		Matrix& apply(const Matrix& m, F f);   //���ض�Ԫ��ִ��f(x,m)

		Matrix& operator=(const T& value); //������ֵ

		Matrix& operator+=(const T& value); //������
		Matrix& operator-=(const T& value); 
		Matrix& operator*=(const T& value); 
		Matrix& operator/=(const T& value); 
		Matrix& operator%=(const T& value); 

	
		Matrix& operator+=(const Matrix& m);//�����
	
		Matrix& operator-=(const Matrix& m);//�����
		
	private:
		Matrix_slice<N> desc;   //����N��Ԫ�ش�С����Ƭ
		std::vector<T> elems;        //Ԫ��
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


