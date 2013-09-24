#include <cmath>

namespace lawa {

template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(XType /*xtype*/, int j, long /*k*/) const
{
    return 1./std::sqrt(1+pow2i<T>(2*j));
}

template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(const Index1D &index) const
{
    return this->operator()(index.xtype,index.j,index.k);
}

}   // namespace lawa

