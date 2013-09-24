#include <bench/aux/random.h>

#include <cstdlib>
#include <ctime>

void
randomInit()
{
    srand((unsigned)(time(0)));
}

template <typename T>
T
randomValue()
{
    return rand()/(T(RAND_MAX)+1);
}
