#ifndef THREADPOOL_HPP
#define THREADPOOL_HPP
#include<iostream>
#include<cstdint>
#include<future>
#include<vector>
#include<queue>
#include<thread>
#include<mutex>
#include<condition_variable>
#include<atomic>
#include<ctime>
using namespace std;
using namespace std::chrono_literals;

template<typename x_value_t,typename y_value_t,typename z_value_t>
struct state_t {
    x_value_t x;
    y_value_t y;
    z_value_t z;
};

template<typename R,typename S,typename T>
void status() {
    typedef atomic<state_t<R, S, T>> atomic_state_t;
    cout << sizeof(atomic_state_t) << "\t" << atomic_state_t().is_lock_free() << endl;
}


class ThreadPool {
private:
    vector<thread> threads;
    queue<function<void(void)>> tasks;

    std::mutex mutex;
    condition_variable cv,cv_wait;

    bool stop_pool;
    //uint32_t active_threads;
    atomic<uint32_t> active_threads;
    const uint32_t capacity;

    template<typename Func, typename...Args,
        typename Rtrn = typename result_of<Func(Args...)>::type>
        auto make_task(Func&& func, Args&&...args)->packaged_task<Rtrn(void)> {
        auto aux = bind(forward<Func>(func), forward<Args>(args)...);//without argument returning func(arg0,arg1...)=aux(void)
        auto task = packaged_task<Rtrn(void)>(aux);
        return task;
    }

    void before_task_hook() {
        active_threads++;
        if (active_threads == 0 && tasks.empty())
        {
            stop_pool = true;
            cv_wait.notify_one();
        }
    }

    void after_task_hook() {
        active_threads--;
    }

public:
    ThreadPool(uint32_t capacity_);

    ~ThreadPool();

    template<typename Func, typename...Args,
        typename Rtrn = typename result_of<Func(Args...)>::type>
        auto enqueue(Func&& func, Args&&...args)->future<Rtrn> {
        auto task = make_task(func, args...);
        auto future = task.get_future();
        auto task_ptr = make_shared<decltype(task)>(move(task));
        {
            lock_guard<std::mutex> lock_guard(mutex);
            if (stop_pool)
                throw runtime_error("enqueue no stopped ThreadPool");

            auto payload = [task_ptr]()->void {
                task_ptr->operator()();
            };
            tasks.emplace(payload);
        }
        cv.notify_all();
        return future;
    }

    void wait_and_stop();

    template<typename Func,typename ...Args>
    void spawn(Func && func,Args&& ...args) {
        if (active_threads < capacity)
            enqueue(func, args...);
        else
            func(args...);
    }
};





void test1();

void test2();

void test3();

#endif 