namespace lawa {

template <typename Basis1D>
InitialCondition1D<Basis1D>::InitialCondition1D(const Function<T> &_payoffInitialCondition,
                                                const Basis1D &_basis, T _left, T _right)
    : payoffInitialCondition(_payoffInitialCondition), basis(_basis), left(_left), right(_right),
      integral(payoffInitialCondition,basis,left,right)
{
    integral.quadrature.setOrder(20);
}

template <typename Basis1D>
typename Basis1D::T
InitialCondition1D<Basis1D>::operator()(XType _e1, int _j1, long _k1, int _deriv1) const
{
    return integral(_j1,_k1,_e1,_deriv1);
}



}   // namespace lawa
