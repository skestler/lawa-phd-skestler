namespace lawa {

template <typename T, typename Basis1D>
FinanceOperator1D<T, BlackScholes, Basis1D>::FinanceOperator1D
                                     (const Basis1D& _basis,
                                      const ProcessParameters1D<T,BlackScholes> &_processparameters,
                                      T _R1, T _R2, int /*order*/,
                                      const int /*internal_compression_level*/)
    : basis(_basis), processparameters(_processparameters),
      R1(_R1), R2(_R2),
      OneDivR2pR1(1./(R2+R1)), OneDivR2pR1squared(OneDivR2pR1*OneDivR2pR1),
      integral(basis,basis)

{

}

template <typename T, typename Basis1D>
T
FinanceOperator1D<T, BlackScholes, Basis1D>::operator()(XType xtype1, int j1, int k1,
                                                        XType xtype2, int j2, int k2) const
{
    T sigma = processparameters.sigma;

    return 0.5*sigma*sigma*( OneDivR2pR1*OneDivR2pR1 * integral(j1,k1,xtype1,1,j2,k2,xtype2,1)
                                       +OneDivR2pR1 * integral(j1,k1,xtype1,0,j2,k2,xtype2,1) );

}

}   // namespace lawa
