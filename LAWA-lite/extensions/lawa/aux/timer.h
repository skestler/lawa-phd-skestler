#ifndef LAWA_AUX_TIMER_H
#define LAWA_AUX_TIMER_H 1

#include <ctime>

namespace lawa {
    
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

} // namespace lawa

#endif //LAWA_AUX_TIMER_H

