#include <cmath>
#include <iostream>
#include <limits>
#include <bench/aux/random.h>

template <typename V>
Axpy<V>::Axpy(long _N)
    : N(_N), x(N), y(N)
{
    for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        x(i) = randomValue<ElementType>();
    }
    for (IndexType i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        y(i) = randomValue<ElementType>();
    }
    alpha = randomValue<ElementType>();
}

template <typename V>
void
Axpy<V>::run(long numberOfComputations)
{
    for (long i=0; i<numberOfComputations; ++i) {
        compute();
    }
}

template <typename V>
long
Axpy<V>::numBaseOperations()
{
    return 2*N;
}

template <typename V>
bool
Axpy<V>::test()
{
#ifdef HAVE_CBLAS
    Vector          xRef = x, yRef = y;
    ElementType     alphaRef = alpha;
    ElementType     tol = std::numeric_limits<ElementType>::epsilon();

    cxxblas::axpy(yRef.length(), alphaRef,
                  xRef.engine().data(), xRef.engine().stride(),
                  yRef.engine().data(), yRef.engine().stride());

    compute();

    for (IndexType i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        if (abs(y(i)-yRef(i))>tol) {
            return false;
        }
    }
    return true;
#else 
    std::cerr << "link with native BLAS for testing" << std::endl;
    return false;
#endif
}
