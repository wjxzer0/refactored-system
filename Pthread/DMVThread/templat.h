#pragma once
#include<iostream>
#include<thread>
#include<vector>
#include<ctime>
#include<math.h>
#include<mutex>
#include<chrono>
#include<condition_variable>
#include<future>
using namespace std;
using namespace chrono_literals;

template<typename value_t, typename index_t>
void init(vector<value_t>& A, vector<value_t>& x, index_t m, index_t n) {
	for (index_t row = 0; row < m; row++)
		for (index_t col = 0; col < n; col++)
			A[row * n + col] = row >= col ? 1 : 0;
	for (index_t col = 0; col < m; col++)
		x[col] = col + 1;
}

template<typename value_t, typename index_t>
void cyclic_parallel_mult(vector<value_t>& A, vector<value_t>& x, vector<value_t>& b,
	index_t m, index_t n, index_t num_threads = 8) {
	auto cyclic = [&](const index_t& id)->void {
		for (index_t row = id; row < m; row += num_threads) {
			value_t accum = value_t(0);
			for (index_t col = 0; col < n; col++)
				accum += A[row * n + col] * x[col];
			b[row] = accum;
		}
	};
	vector<thread> threads;
	for (index_t id = 0; id < num_threads; id++)
		threads.emplace_back(cyclic, id);
	for (auto& thread : threads)
		thread.join();

}

template<typename value_t, typename index_t>
void sequential_mult(vector<value_t>& A, vector<value_t>& x, vector<value_t>& b,
	index_t m, index_t n) {
	for (index_t row = 0; row < m; row ++) {
		value_t accum = value_t(0);
		for (index_t col = 0; col < n; col++)
			accum += A[row * n + col] * x[col];
		b[row] = accum;
	}
}

template<typename value_t, typename index_t>
void block_cyclic_parallel_mult(vector<value_t>& A, vector<value_t>& x, vector<value_t>& b,
	index_t m, index_t n, index_t num_threads = 8, index_t chunk_size = 64 / sizeof(value_t)) {
	auto block_cyclic = [&](const index_t& id)->void {
		const index_t offset = id * chunk_size;
		const index_t stride = num_threads * chunk_size;

		for (index_t lower = offset; lower < m; lower += stride) {
			const index_t upper = min(lower + chunk_size, m);

			for (index_t row = lower; row < upper; row++) {
				value_t accum = value_t(0);
				for (index_t col = 0; col < n; col++)
					accum += A[row * n + col] * x[col];
				b[row] = accum;
			}
		}
	};
	vector<thread> threads;
	for (index_t id = 0; id < num_threads; id++)
		threads.emplace_back(block_cyclic, id);
	for (auto& thread : threads)
		thread.join();

}

template<typename value_t, typename index_t>
void printb(vector<value_t>& b, index_t m) {
	for (index_t row = 0; row < m; row++) {
		cout << b[row] << endl;
	}
}

struct pack_t {
	uint64_t ying;
	uint64_t yang;

	pack_t() :ying{ 0 }, yang{ 0 }{}
};

void sequential_increment(volatile pack_t& pack);

void false_sharing_increment(volatile pack_t& pack);

void test1();

void test2();

void test3();//sleep student

void test4();//pong ,ping

void test5();

void test6();