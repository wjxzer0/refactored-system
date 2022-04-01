#include"threadpool.hpp"


void ThreadPool::wait_and_stop() {
    unique_lock<std::mutex> unique_lock(mutex);
    auto predicate = [&]()->bool{
        return stop_pool;
    };
    cv_wait.wait(unique_lock, predicate);
}

void test1() {
    cout << "in" << endl;
    //this_thread::sleep_for(3s);
	ThreadPool PT(8);

    cout << "go" << endl;
	auto squar = [](const uint64_t x) {return x * x; };
	const uint64_t num_task = 32;
	vector<future<uint64_t>> futures;
	for (uint64_t task = 0; task < num_task; task++) {
		auto future = PT.enqueue(squar, task);
		futures.emplace_back(move(future));
	}
	for (auto& future : futures) {
		cout << future.get() << endl;
	}
}


ThreadPool::~ThreadPool() {
	{
		lock_guard<std::mutex> lock_guard(mutex);
		stop_pool = true;
	}

	cv.notify_all();
	for (auto& thread : threads)
		thread.join();
}

ThreadPool::ThreadPool(uint32_t capacity_) :stop_pool{ false }, active_threads{ 0 }, capacity{ capacity_ }{
    auto wait_loop = [this]()->void {
        while (true) {
            function<void(void)>task;

            {
                std::mutex mutextmp;
                unique_lock<std::mutex> unique_lock(mutextmp);
                //unique_lock<std::mutex> unique_lock(mutex);
                cout << "in!" << endl;
                auto predicate = [this]()->bool {return (stop_pool) || !(tasks.empty()); };
                cv.wait(unique_lock, predicate);
                if (stop_pool && tasks.empty())
                    return;

                task = move(tasks.front());
                tasks.pop();
                before_task_hook();
            }

            task();
            {
                lock_guard<std::mutex> lock_guard(mutex);
                after_task_hook();
            }

        }
    };

    for (uint64_t id = 0; id < capacity; id++)
        threads.emplace_back(wait_loop);
};

void test2() {
    mutex mutex;
    vector<thread> threads;
    const uint64_t num_threads = 10;
    const uint64_t num_iters = 100'000'000;
    /*auto lock_count = [&](volatile uint64_t* counter, const auto& id)->void {
        for (uint64_t i = id; i < num_iters; i += num_threads) {
            lock_guard<std::mutex> lock_guard(mutex);
            (*counter)++;
        }
    };*/
    auto atomic_count = [&](volatile atomic<uint64_t>* counter, const auto& id)->void {
        for (uint64_t i = id; i < num_iters; i += num_threads) {
            auto previous = counter->load();
            while (previous < i && !counter->compare_exchange_weak(previous, i)) {}
        }
    };
    clock_t start;
    clock_t end;
    
    atomic<uint64_t> atomic_counter{ 0 };
    threads.clear();
    start = clock();
    for (uint64_t id = 0; id < num_threads; id++)
        threads.emplace_back(atomic_count, &atomic_counter, id);
    for (auto& thread : threads)
        thread.join();
    end = clock();
    cout << end - start << "time" << endl;


    /*uint64_t counter = 0;
    threads.clear();
    start = clock();
    for (uint64_t id = 0; id < num_threads; id++)
        threads.emplace_back(lock_count, &counter, id);
    for (auto& thread : threads)
        thread.join();
    end = clock();
    cout << end - start <<"time" << endl;*/
    
    

    cout <<  atomic_counter << endl;

}

void test3() {
    cout << "size\tlock_free?" << endl;

    status<uint8_t, uint8_t, uint8_t>();
    status<uint16_t, uint8_t, uint8_t>(); 
    status<uint16_t, uint16_t, uint8_t>();
    status<uint32_t, uint16_t, uint16_t>();
    status<uint32_t, uint32_t, uint16_t>();
    status<uint64_t, uint32_t, uint32_t>();

}