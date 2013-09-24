#include <cassert>

namespace lawa {

template <typename T, FunctionSide Side, Construction Cons>
T
BasisFunction<T,Side,Periodic,Cons>::operator()(T /*x*/, int /*j*/, long /*k*/, 
                                                unsigned short /*deriv*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, Construction Cons>
PeriodicSupport<T>
BasisFunction<T,Side,Periodic,Cons>::support(int /*j*/, long /*k*/) const 
{
    assert(0);
    return PeriodicSupport<T>();
}

template <typename T, FunctionSide Side, Construction Cons>
DenseVector<Array<T> >
BasisFunction<T,Side,Periodic,Cons>::singularSupport(int /*j*/, long /*k*/) const 
{
    assert(0);
    return DenseVector<Array<T> >(); 
}

template <typename T, FunctionSide Side, Construction Cons>
T
BasisFunction<T,Side,Periodic,Cons>::tic(int /*j*/) const
{
    assert(0);
    return 0.;
}

} // namespace lawa

