#pragma once
namespace M {
template<size_t N>
struct Matrix_slice {
	Matrix_slice() = default;       //空矩阵：无元素

	Matrix_slice(size_t s, std::initializer_list<size_t> exts);
	Matrix_slice(size_t s, std::initializer_list<size_t> exts, std::initializer_list<size_t> strs);

	template<typename... Dims>
	Matrix_slice(Dims... dims);    //N个维度大小

	template<typename... Dims,typename=std::enable_if<All(is_convertible<Dims,size_t>()...)>>
	size_t operator()(Dims... dims)const {
		static_assert(sizeof...(Dims) == N, "");
		size_t args[N]{ size_t(dims...) };
		return 0;
		//inner_product(args, args + N, strides.begin(), size_t(0));
	}    //从一组下标计算索引
	
	size_t size;  //元素总数
	size_t start;  //起始偏移量
	size_t* extents[N];//array<size_t, N> extents;  //每个维度大小
	size_t* strides[N];//array<size_t, N> strides;   // 每个维度上元素间的偏移量
};
}