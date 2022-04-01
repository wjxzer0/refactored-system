#pragma once
#include<iostream>
#include"Matrix_slice.h"
#include"slice.h"
using namespace std;
using namespace M;
namespace Matrix_impl {

	template<size_t N, typename l, typename List>
	enable_if<(N > 1), void>add_extends(l& first, const List& list) {
		assert(check_non_jagged(list));
		*first = list.size();
		add_extends<N - 1>(++first, *list.begin());
	}

	template<size_t N, typename l, typename List>
	enable_if<(N == 1), void>add_extends(l& first, const List& list) {
		*first++ = list.size();
	}


	template<size_t N, typename List>
	array<size_t, N> derive_extents(const List& list) {
		array<size_t, N> a;
		auto f = a.begin();
		add_extends<N>(f, list);       //将维度大小添加到a中
		return a;
	}

	template<typename List>
	bool check_non_jagged(const List& list) {//检查每行是否包含相同数量的元素
		auto i = list.begin();
		for (auto j = i + 1; j != list.end(); ++j)
			if (i->size() != j->size())
				return false;
		return true;
	}


	template<typename T,typename Vec>
	void insert_flat(initializer_list<T>list, Vec& vec) {
		add_list(list.begin(), list.end(), vec);
	}


	//嵌套的initializer_list初始化
	template<typename T,typename Vec>//接受initializer_list
	void add_list(const initializer_list<T>* first, const initializer_list<T>* last, Vec& vec) {
		for (; first != last; ++first)
			add_list(first->begin(), first->end(), vec);
	}

	template<typename T, typename Vec>//元素插入
	void add_list(const T* first, const T* last, Vec& vec) {
		vec.insert(vec.end(), first, last);
	}

	template<typename...Args>
	constexpr bool Requesting_element() {
		return All(is_convertible<Args, size_t>()...);
	}

	constexpr bool All() { return true; }

	template<typename...Args>
	constexpr bool All(bool b, Args... args) {
		return b && All(args...);
	}

	constexpr bool Some() { return true; }

	template<typename...Args>
	constexpr bool Some(bool b, Args... args) {
		return b || All(args...);
	}

	template<typename...Args>
	constexpr bool Requesting_slice() {
		return All((is_convertible<Args, size_t>() || is_same<Args, slice>())...)
			&& Some(is_same<Args,slice>()...);
	}

	template<size_t N,typename...Dims>
	bool check_bounds(const Matrix_slice<N>& slice, Dims...dims) {
		size_t indexes[N]{ size_t(dims)... };
		return equal(indexes, indexes + N, slice.extends, less<size_t>{});
	}

    //initializer
	template<typename T, size_t N>
	struct Matrix_init {
		using type = initializer_list<typename Matrix_init<T, N - 1>::type>;
	};

	template<typename T>
	struct Matrix_init<T, 1> {
		using type = initializer_list<T>;
	};

	template<typename T>
	struct Matrix_init<T, 0>;

	template<size_t N,typename T,typename...Args>
	size_t do_slice(const Matrix_slice<N>& os, Matrix_slice<N>& ns, const T& s, const Args&...args) {
		size_t m = do_slice_dim(os, ns, s);                 //do_slice_dim<sizeof...(Args) + 1>(os, ns, s);
		size_t n = do_slice(os, ns, args...);
		return m + n;
	}

	template<size_t N>
	size_t do_slice(const Matrix_slice<N>& os, const Matrix_slice<N>& ns) {
		return 0;
	}

	template<size_t N,typename T, typename...Args>
	size_t do_slice_dim(const Matrix_slice<N>& os, Matrix_slice<N>& ns, const T& s) {
		//...
		return 1;
	}

	template<size_t N,typename T>
	size_t slice_dim(size_t n, const Matrix_slice<N> desc, const Matrix_slice<N>& row) {
		//...
		return 1;
	}
}