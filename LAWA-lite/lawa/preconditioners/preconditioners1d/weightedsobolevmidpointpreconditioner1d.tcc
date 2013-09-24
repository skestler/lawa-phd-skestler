#include <cmath>

namespace lawa {

template <typename T, typename Basis>
WeightedSobolevMidPointPreconditioner1D<T,Basis>::WeightedSobolevMidPointPreconditioner1D
                                                 (const Basis &basis, const Function<T> &weight,
                                                  const int sobolev_order)
    : _basis(basis), _weight(weight), _sobolev_order(sobolev_order),
      _integral(basis,basis)
{
    assert(sobolev_order<=1);
}

template <typename T, typename Basis>
T
WeightedSobolevMidPointPreconditioner1D<T,Basis>::operator()(XType xtype, int j, long k) const
{
    T center = 0.5*(_basis.generator(xtype).support(j,k).l2+_basis.generator(xtype).support(j,k).l1);
    if (_sobolev_order==0) {
        return 1./std::sqrt(_weight(center) * _integral(j,k,xtype,0,j,k,xtype,0));
    }
    else if (_sobolev_order==1) {
        return 1./std::sqrt( _weight(center) * ( _integral(j,k,xtype,0,j,k,xtype,0)
                            + _integral(j,k,xtype,1,j,k,xtype,1) ) );
    }
    else {
        assert(0);
        return 0.;
    }

}

template <typename T, typename Basis>
T
WeightedSobolevMidPointPreconditioner1D<T,Basis>::operator()(const Index1D &index) const
{
    return this->operator()(index.xtype,index.j,index.k);
}

}   // namespace lawa
