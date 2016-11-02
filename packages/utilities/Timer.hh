#ifndef Timer_hh
#define Timer_hh

#include <chrono>

class Timer
{
private:

    bool going_;
    bool timed_;

    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point end_time_;
    std::chrono::duration<double> elapsed_time_;
    
public:

    Timer();

    void start();
    void stop();
    double time();
    void print_time();
};

#endif
