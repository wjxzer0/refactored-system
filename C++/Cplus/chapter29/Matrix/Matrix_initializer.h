#pragma once
namespace M {
	
	template<typename T, size_t N>
	using Matrix_initializer = typename Matrix_impl::Matrix_init<T, N>::type;
}