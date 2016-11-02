#include "Timer.hh"

#include <chrono>
#include <iomanip>
#include <iostream>

#include "Check.hh"

using namespace std;
using namespace std::chrono;

Timer::
Timer()
{
    going_ = false;
    timed_ = false;
}

void Timer::
start()
{
    if (going_)
    {
        CheckMsg(false, "timer already running");
    }
    
    start_time_ = high_resolution_clock::now();
    going_ = true;
}

void Timer::
stop()
{
    if (!going_)
    {
        CheckMsg(false, "timer is not running");
    }

    end_time_ = high_resolution_clock::now();
    going_ = false;
    timed_ = true;

    elapsed_time_ = end_time_ - start_time_;
}

double Timer::
time()
{
    if (!timed_)
    {
        CheckMsg(false, "timer has not run");
    }
    
    return elapsed_time_.count();
}

void Timer::
print_time()
{
    if (!timed_)
    {
        CheckMsg(false, "timer has not run");
    }

    cout << setprecision(10) << elapsed_time_.count() << endl;
}
