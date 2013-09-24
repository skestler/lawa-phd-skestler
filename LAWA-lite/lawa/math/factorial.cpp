#include <lawa/math/factorial.h>

namespace lawa {

long
factorial(int n)
{
    long fac = 1;
    for (int i=2; i<=n; ++i) {
        fac *= i;
    }
    return fac;
}

} // namespace lawa

