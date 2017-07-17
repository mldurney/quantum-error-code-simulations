#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <vector>
#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

class ThreadPool {
public:
	explicit ThreadPool(unsigned int numThreads = 4);
	~ThreadPool();
	void addTask(std::function <void(void)> f);
	void joinThreads(bool waitForAll = true);
	int size() const { return getMaxThreads(); }
	int getMaxThreads() const { return maxThreads; }
	int getTasksRemaining() const { return tasksRemaining; }

private:
	void runTask();
	std::function<void(void)> nextTask();
	void waitAll();

	std::vector<std::thread> threads;
	std::queue<std::function <void(void)>> taskQueue;
	std::atomic_uint maxThreads;
	std::atomic_uint tasksRemaining;
	std::atomic_bool exitThreading;
	std::atomic_bool finished;
	std::condition_variable task_available_var;
	std::condition_variable wait_var;
	std::mutex wait_mutex;
	std::mutex queue_mutex;
};

#endif /* THREADPOOL_H */