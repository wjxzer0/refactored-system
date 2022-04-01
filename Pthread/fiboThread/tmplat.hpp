#pragma once
#include<iostream>
#include<future>
#include<cstdint>
#include <thread>

using namespace std;

template<typename Func, typename...Args,
    typename Rtrn = typename result_of<Func(Args...)>::type>
    auto make_task(Func&& func, Args&&...args)->packaged_task<Rtrn(void)> {
    auto aux = bind(forward<Func>(func), forward<Args>(args)...);//without argument returning func(arg0,arg1...)=aux(void)
    auto task = packaged_task<Rtrn(void)>(aux);
    return task;
}
    
namespace TEST{
    template <typename value_t, typename index_t>
    void fibo(value_t n, promise<uint64_t>&& result)
    {
        value_t a_0 = 0;
        value_t a_1 = 1;
        for (index_t index = 0; index < n; index++)
        {
            const value_t tmp = a_0;
            a_0 = a_1;
            a_1 += tmp;
        }
        result.set_value(a_0);
    }

    void test1();
    void test2();
    void test3();
    uint64_t fibo(uint64_t n);
}



