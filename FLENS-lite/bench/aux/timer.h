#ifndef BENCH_TIMER_H
#define BENCH_TIMER_H 1

#include <ctime>

struct Timer
{
    void
    start();

    void
    stop();

    double
    elapsed();

    clock_t _start;
    clock_t _stop;
};

#endif // BENCH_TIMER_H