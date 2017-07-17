#include "threadpool.h"

ThreadPool::ThreadPool(unsigned int numThreads) : maxThreads(numThreads), exitThreading(false), finished(false) {
	for (unsigned int i = 0; i < maxThreads; ++i) {
		threads.push_back(std::thread([this] {this->runTask(); }));
	}
}

ThreadPool::~ThreadPool() {
	joinThreads();
}

void ThreadPool::runTask() {
	while (!exitThreading) {
		nextTask()();
		--tasksRemaining;
		wait_var.notify_one();
	}
}

std::function<void(void)> ThreadPool::nextTask() {
	std::function<void(void)> task;
	std::unique_lock<std::mutex> lock(queue_mutex);

	task_available_var.wait(lock, 
		[this]()-> bool {return taskQueue.size() || exitThreading; });

	if (!exitThreading) {
		task = taskQueue.front();
		taskQueue.pop();
	}
	else {
		task = [] {};
		++tasksRemaining;
	}

	return task;
}

void ThreadPool::addTask(std::function <void(void)> f) {
	std::lock_guard<std::mutex> guard(queue_mutex);
	taskQueue.push(f);
	++tasksRemaining;
	task_available_var.notify_one();
}

void ThreadPool::joinThreads(bool waitForAll) {
	if (!finished) {
		if (waitForAll) {
			waitAll();
		}

		exitThreading = true;
		task_available_var.notify_all();

		for (auto &t : threads) {
			if (t.joinable()) {
				t.join();
			}
		}

		finished = true;
	}
}

void ThreadPool::waitAll() {
	if (tasksRemaining > 0) {
		std::unique_lock<std::mutex> lock(wait_mutex);
		wait_var.wait(lock, [this] {return this->tasksRemaining == 0; });
		lock.unlock();
	}
}

