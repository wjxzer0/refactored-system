#include"Matrix.h"

namespace M {
	template<typename T, size_t N>
	template<typename... Exis>
	Matrix<T, N>::Matrix(Exis...exts) :desc{ exts... }//����ά�ȴ�С
		, elems(desc.size)//����desc.size��Ԫ�أ��������ǽ���Ĭ�ϳ�ʼ��
	{}



	template<typename T, size_t N>
	Matrix<T, N>::Matrix(Matrix_initializer<T, N> init) {
		Matrix_impl::derive_extents(init, desc.extents); //�ӳ�ʼ�����ж�ά�ȴ�С
		elems.reserve(desc.size());        //Ϊ��Ƭ�����ռ�
		Matrix_impl::insert_flat(init, elems);    //�ó�ʼ�����б��ʼ��
		assert(elems.size() == desc.size);
	}

	template<typename T, size_t N>
	template<typename U>
	Matrix<T, N>::Matrix(const Matrix_ref<U, N>& x)
		:desc{ x.desc }, elems{ x.begin(),x.end() }//����descԪ��
	{
		static_assert(is_convertible<U, T>(), "Matrix constructor:incompatible element types");
	}

	template<typename T, size_t N>
	template<typename U>
	Matrix<T, N>& Matrix<T, N>::operator=(const Matrix_ref<U, N>& x) {
		static_assert(is_convertible<U, T>(), "Matrix constructor:incompatible element types");

		desc = x.desc;
		elems.assign(x.begin(), x.end());
		return *this;
	}

	//��������
	template<typename T, size_t N>
	template<typename F>
	Matrix<T, N>& Matrix<T, N>::apply(F f) {
		for (auto& x : elems)f(x);
		return *this;
	}

	template<typename T, size_t N>
	Matrix<T, N>& Matrix<T, N>::operator+=(const T& value) {
		return apply([&](T& a) {a += value; });
	}

	template<typename T, size_t N>
	Matrix<T, N>& Matrix<T, N>::operator-=(const T& value) {
		return apply([&](T& a) {a -= value; });
	}

	template<typename T, size_t N>
	Matrix<T, N>& Matrix<T, N>::operator*=(const T& value) {
		return apply([&](T& a) {a *= value; });
	}

	template<typename T, size_t N>
	Matrix<T, N>& Matrix<T, N>::operator/=(const T& value) {
		return apply([&](T& a) {a /= value; });
	}

	template<typename T, size_t N>
	Matrix<T, N>& Matrix<T, N>::operator%=(const T& value) {
		return apply([&](T& a) {a %= value; });
	}

	//�������

	template<typename T, size_t N>
	 Matrix<T, N>& Matrix<T, N>::operator+=(const Matrix& m) {
		static_assert(m.order() == N, "+=:mismatched Matrix dimentions");
		assert(same_extents(desc, m.descriptor()));

		return apply(m, [](T& a, Matrix<T,N>& b) {a += b; });
	}

	 template<typename T, size_t N>
	 Matrix<T, N>& Matrix<T, N>::operator-=(const Matrix& m) {
		 static_assert(m.get_order() == N, "+=:mismatched Matrix dimentions");
		 assert(same_extents(desc, m.descriptor()));

		 return apply(m, [](T& a, Matrix<T, N>& b) {a -= b; });
	 }

	template<typename T, size_t N>
	template<typename F>
	Matrix<T, N>& Matrix<T, N>::apply(const Matrix<T,N>& m, F f) {
		assert(same_extents(desc, m.descriptor()));
		for (auto i = this.begin(), j = m.end(); i != this.end(); ++i, ++j)
			f(*i, *j);
		return *this;
	}

	template<typename T,size_t N>
	Matrix_ref<T, N - 1>Matrix<T, N>::row(size_t n) {
		//assert(n < row());
		Matrix_slice<N - 1> row;
		Matrix_impl::slice_dim<0>(n, desc, row);
		return { row,data() };
	}

	template<typename T, size_t N>
	Matrix_ref<const T, N - 1>Matrix<T, N>::row(size_t n)const {
		//assert(n < row());
		Matrix_slice<N - 1> row;
		Matrix_impl::slice_dim<0>(n, desc, row);
		return { row,data() };
	}

	template<typename T, size_t N>
	Matrix_ref<T, N - 1>Matrix<T, N>::col(size_t n) {
		//assert(n < col());
		Matrix_slice<N - 1> col;
		Matrix_impl::slice_dim<1>(n, desc, col);
		return { col,data() };
	}

	template<typename T, size_t N>
	Matrix_ref<const T, N - 1>Matrix<T, N>::col(size_t n)const {
		//assert(n < row());
		Matrix_slice<N - 1> col;
		Matrix_impl::slice_dim<1>(n, desc, col);
		return { col,data() };
	}

	template<typename T, size_t N>
	template<typename...Args>
	enable_if<Matrix_impl::Requesting_element<Args...>(), T&>
		Matrix<T, N>::operator()(Args... args) {
		assert(Matrix_impl::check_bounds(desc, args...));
		return *(data() + desc(args...));
	}

	template<typename T, size_t N>
	template<typename...Args>
	enable_if<Matrix_impl::Requesting_element<Args...>(),const T&>
		Matrix<T, N>::operator()(Args... args)const {
		assert(Matrix_impl::check_bounds(desc, args...));
		return *(data() + desc(args...));
	}

	template<typename T, size_t N>
	template<typename...Args>
	enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T, N>>
		Matrix<T, N>::operator()(const Args... args) {
		Matrix_slice<N>d;
		d.start = Matrix_impl::do_slice(desc, d, args...);
			return { d,data() };
	}

	template<typename T, size_t N>
	template<typename...Args>
	enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T, N>>
		Matrix<T, N>::operator()(const Args... args)const {
		Matrix_slice<N>d;
		d.start = Matrix_impl::do_slice(desc, d, args...);
			return { d,data() };
	}

}



