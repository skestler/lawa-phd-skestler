#include <cmath>

namespace lawa {

template <typename T, typename Basis>
WeightedSobolevNormPreconditioner1D<T,Basis>::WeightedSobolevNormPreconditioner1D
                                              (const Basis &basis, const Function<T> &weight,
                                               const int sobolev_order)
    : _sobolev_order(sobolev_order), _integral(weight,basis,basis)
{
    assert(sobolev_order<=1);
    _integral.quadrature.setOrder(8*basis.d);
}

template <typename T, typename Basis>
T
WeightedSobolevNormPreconditioner1D<T,Basis>::operator()(XType xtype, int j, long k) const
{
    if (_sobolev_order==0) {
        return 1./std::sqrt(_integral(j,k,xtype,0,j,k,xtype,0));
    }
    else if (_sobolev_order==1) {
        return 1./std::sqrt(  _integral(j,k,xtype,0,j,k,xtype,0)
                            + _integral(j,k,xtype,1,j,k,xtype,1));
    }
    else {
        assert(0);
        return 0.;
    }

}

template <typename T, typename Basis>
T
WeightedSobolevNormPreconditioner1D<T,Basis>::operator()(const Index1D &index) const
{
    return this->operator()(index.xtype,index.j,index.k);
}

}   // namespace lawa

