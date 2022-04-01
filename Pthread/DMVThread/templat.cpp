#include"templat.h"

void test1() {
	clock_t start;
	clock_t end;
	start = clock();
	const uint64_t m = 10000;
	const uint64_t n = 10000;


	vector<uint64_t>A(m * n);
	vector<uint64_t>x(n);
	vector<uint64_t>b(m);
	end = clock();
	cout << end - start << endl;


	init(A, x, m, n);

	start = clock();

	cout << start - end << endl;

	cyclic_parallel_mult(A, x, b, m, n);

	end = clock();
	cout <<"cyclic_parallel_mult :  " << end - start << endl;
	
	start = clock();
	sequential_mult(A, x, b, m, n);
	end = clock();
	cout << "sequrntial_mult :   " << end - start << endl;

	start = clock();
	block_cyclic_parallel_mult(A, x, b, m, n);
	end = clock();
	cout << "sequrntial_mult :   " << end - start << endl;

}

void sequential_increment(volatile pack_t& pack) {
	for (uint64_t index = 0; index < 1UL << 20; index++) {
		pack.ying++;
		pack.yang++;
	}
}

void false_sharing_increment(volatile pack_t& pack) {
	auto eval_ying = [&pack]()->void {
		for (uint64_t index = 0; index < 1UL << 20; index++) 
			pack.ying++;
	};

	auto eval_yang = [&pack]()->void {
		for (uint64_t index = 0; index < 1UL << 20; index++)
			pack.yang++;
	};

	thread ying_thread(eval_ying);
	thread yang_thread(eval_yang);
	ying_thread.join();
	yang_thread.join();
}

void test2() {
	pack_t seq_pack;
	clock_t start;
	clock_t end;

	start = clock();
	sequential_increment(seq_pack);
	end = clock();
	cout << end - start << endl;

	cout << seq_pack.ying << "  " << seq_pack.ying << endl;

	pack_t par_pack;

	start = clock();
	false_sharing_increment(par_pack);
	end = clock();
	cout << end - start << endl;

	cout << par_pack.ying << "  " << par_pack.ying << endl;

}

void test3() {
	mutex mutex;
	condition_variable cv;
	bool time_for_breakfast = false;

	auto student = [&]()->void {
		{
			unique_lock<std::mutex> unique_lock(mutex);//lock

			while (!time_for_breakfast)//lock is released during wait
				cv.wait(unique_lock);//cv.wait(unique_lock,[&](){return time_for_breakfast;});
		};//released lock

		cout << "Time to make some coffee!" << endl; 
	};

	thread myshread(student);
	this_thread::sleep_for(2s);

	{
		lock_guard<std::mutex>lock_guard(mutex);
		time_for_breakfast = true;
	}
	//time_for_breakfast = true;
	cout << "ture!" << endl;
	cv.notify_all();
	
	cout << "??" << endl;
	myshread.join();

	
}

void test4() {
	mutex mutex;
	std::mutex mutex2;
	condition_variable cv;
	bool is_ping = false;
	auto ping = [&](int id)->void {
		while (true) {
			cout << "??i" << endl;
			unique_lock<std::mutex> unique_lock(mutex);
			cout << "get lock: " << id << endl;
			
			cv.wait(unique_lock, [&]() {return is_ping; });

			this_thread::sleep_for(2s);
			cout << "ping  " <<id<< endl;

			is_ping = !is_ping;
			//cv.notify_one();
		}
	};

	auto pong = [&]()->void {
		while (true) {
			//cout << "??o" << endl;
			//unique_lock<std::mutex> unique_lock(mutex);
			cout << "pong lock" << endl;
			//cv.wait(unique_lock, [&]() {return !is_ping; });

			
			cout << "pong" << endl;

			
			cv.notify_all();
			this_thread::sleep_for(10s);
			cout << "notify ping" << endl;
			is_ping = !is_ping;
			cout << is_ping << endl;
			//cv.notify_all();

		}
	};
	
	thread ping_thread(ping,1);
	thread pong_thread(pong);
	thread ping_thread2(ping,2);
	this_thread::sleep_for(10s);
	cout << "ping ping pong pong!" << endl;
	ping_thread.join();
	pong_thread.join();
	ping_thread2.join();

}

void test5() {
	promise<void> promise;
	auto future = promise.get_future().share();
	auto student = [&]()->void {
		future.get();
		cout << "time to breakfast! " << endl;
	};
	thread mythread(student);
	thread mythread2(student);
	this_thread::sleep_for(2s);
	promise.set_value();

	mythread.join();
	mythread2.join();
}