#include<cstdint>
#include <thread>
#include<vector>
#include"tmplat.hpp"

void TEST::test1() {
    const uint64_t num_threads = 32;
    vector<thread> threads;

    vector<future<uint64_t>> results;

    for (uint64_t id = 0; id < num_threads; id++)
    {
        promise<uint64_t> promise;
        results.emplace_back(promise.get_future());
        threads.emplace_back(fibo<uint64_t, uint64_t>, id, move(promise));
    }

    for (auto& result : results)
    {
        cout << result.get() << endl;
    }

    for (auto& thread : threads)
        thread.join();
}

uint64_t TEST::fibo(uint64_t n)
{
    uint64_t a_0 = 0;
    uint64_t a_1 = 1;
    for (uint64_t index = 0; index < n; index++)
    {
        const uint64_t tmp = a_0;
        a_0 = a_1;
        a_1 += tmp;
    }
    return a_0;
}

void TEST::test2() {
    const uint64_t num_threads = 32;
    vector<thread> threads;
    vector<future<uint64_t>> results;
    for (uint64_t id = 0; id < num_threads; id++)
    {
        auto task = make_task<uint64_t (&)(uint64_t)>(fibo, id);
        results.emplace_back(task.get_future());
        threads.emplace_back(std::move(task));
    }
    for (auto& result : results)
    {
        cout << result.get() << endl;
    }
    for (auto& thread : threads)
        thread.join();
}

void TEST::test3() {
    const uint64_t num_threads = 32;
    vector<future<uint64_t>> results;
    for (uint64_t id = 0; id < num_threads; id++)
    {
        results.emplace_back(async<uint64_t(&)(uint64_t)>(launch::async, fibo, id));
    }
    for (auto& result : results)
    {
        cout << result.get() << endl;
    }
}






    
