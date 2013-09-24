#include <lawa/aux/timer.h>

namespace lawa {

void
Timer::start()
{
    _start = clock();
}

void
Timer::stop()
{
    _stop = clock();
}

double
Timer::elapsed()
{
    return double(_stop - _start) / CLOCKS_PER_SEC;
}

} // namespace lawa

