#pragma once

// Adaptation of implementation taken from:
// http://stackoverflow.com/questions/12805041/c-equivalent-to-javas-blockingqueue

#include <condition_variable>
#include <functional>
#include <deque>
#include <mutex>

class queue_closed_exception : std::exception {
public:
    virtual ~queue_closed_exception() {}
    virtual const char *what() const noexcept { return "Cannot push to closed queue."; }
};

template <typename T>
class queue_t
{
private:
    std::mutex              d_mutex;
    std::condition_variable d_condition;
    std::deque<T>           d_queue;
    bool		    d_closed = false;

public:
    void push(T const& value) {
        {
            std::unique_lock<std::mutex> lock(d_mutex);

	    if(d_closed) {
		throw queue_closed_exception();
	    }

            d_queue.push_front(value);
        }
	d_condition.notify_one();
    }

    bool pop(T &result) {
	// Obtain lock
        std::unique_lock<std::mutex> lock(this->d_mutex);

	// Block until queue has item or is closed.
        this->d_condition.wait(lock, [=]{ return !this->d_queue.empty() || d_closed; });

	// Return item or communicate closed
	if(!d_queue.empty()) {
	    result = this->d_queue.back();
	    d_queue.pop_back();

	    return true;
	} else {
	    // queue is closed.
	    return false;
	}
    }

    void close() {
        {
            std::unique_lock<std::mutex> lock(d_mutex);
	    if(d_closed) {
		throw queue_closed_exception();
	    }
	    d_closed = true;
	}
        d_condition.notify_all();
    }

    size_t size() {
	std::unique_lock<std::mutex> lock(d_mutex);
	return d_queue.size();
    }
};
