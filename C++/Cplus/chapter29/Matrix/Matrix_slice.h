#pragma once
namespace M {
template<size_t N>
struct Matrix_slice {
	Matrix_slice() = default;       //�վ�����Ԫ��

	Matrix_slice(size_t s, std::initializer_list<size_t> exts);
	Matrix_slice(size_t s, std::initializer_list<size_t> exts, std::initializer_list<size_t> strs);

	template<typename... Dims>
	Matrix_slice(Dims... dims);    //N��ά�ȴ�С

	template<typename... Dims,typename=std::enable_if<All(is_convertible<Dims,size_t>()...)>>
	size_t operator()(Dims... dims)const {
		static_assert(sizeof...(Dims) == N, "");
		size_t args[N]{ size_t(dims...) };
		return 0;
		//inner_product(args, args + N, strides.begin(), size_t(0));
	}    //��һ���±��������
	
	size_t size;  //Ԫ������
	size_t start;  //��ʼƫ����
	size_t* extents[N];//array<size_t, N> extents;  //ÿ��ά�ȴ�С
	size_t* strides[N];//array<size_t, N> strides;   // ÿ��ά����Ԫ�ؼ��ƫ����
};
}