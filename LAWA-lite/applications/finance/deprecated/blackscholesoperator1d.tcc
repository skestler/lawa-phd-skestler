namespace lawa {

template <typename T, typename Basis1D>
FinanceOperator1D<T, BlackScholes, Basis1D>::FinanceOperator1D
                                     (const Basis1D& _basis,
                                      const ProcessParameters1D<T,BlackScholes> &_processparameters,
                                      const T _eta, T _R1, T _R2, int order,
                                      const int /*internal_compression_level*/)
    : basis(_basis), processparameters(_processparameters),
      eta(_eta), exponentialweightfunction(),
      weight(exponentialweightfunction.weight,exponentialweightfunction.singularPoints),
      dweight(exponentialweightfunction.dweight,exponentialweightfunction.singularPoints),
      R1(_R1), R2(_R2),
      OneDivR2pR1(1./(R2+R1)), OneDivR2pR1squared(OneDivR2pR1*OneDivR2pR1),
      integral(basis,basis), integral_weight(weight, basis, basis, -R1, R2),
      integral_dweight(dweight, basis, basis, -R1, R2)
{
     integral_weight.quadrature.setOrder(order);
     integral_dweight.quadrature.setOrder(order);
     exponentialweightfunction.setEta(eta);
}

template <typename T, typename Basis1D>
T
FinanceOperator1D<T, BlackScholes, Basis1D>::operator()(XType xtype1, int j1, int k1,
                                                        XType xtype2, int j2, int k2) const
{
    T sigma = processparameters.sigma;
    if (eta!=0) {
        return  0.5*sigma*sigma*(
                  integral_weight(j1,k1,xtype1,1,j2,k2,xtype2,1)
                + integral_weight(j1,k1,xtype1,0,j2,k2,xtype2,1)
                + integral_dweight(j1,k1,xtype1,0,j2,k2,xtype2,1) );
    }
    else {
        return 0.5*sigma*sigma*( OneDivR2pR1*OneDivR2pR1 * integral(j1,k1,xtype1,1,j2,k2,xtype2,1)
                                           +OneDivR2pR1 * integral(j1,k1,xtype1,0,j2,k2,xtype2,1) );
    }


}

}   // namespace lawa
