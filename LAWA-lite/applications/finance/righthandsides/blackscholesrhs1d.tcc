namespace lawa {

template <typename T, typename Basis1D>
OptionRHS1D<T, Put, BlackScholes, Basis1D>::OptionRHS1D
                                     (const OptionParameters1D<T,Put> &_optionparameters,
                                      const ProcessParameters1D<T,BlackScholes> &_processparameters,
                                      const Basis1D &_basis, T _R1, T _R2, bool _excessToPayoff)
: optionparameters(_optionparameters), processparameters(_processparameters),
  basis(_basis),
  R1(_R1), R2(_R2),
  OneDivSqrtR2pR1(1./std::sqrt(R2+R1)), OneDivR2pR1(1./(R2+R1)), R1DivR1pR2(R1/(R1+R2)),
  excessToPayoff(_excessToPayoff)
{
    if (R1==0 && R2==1) {
        OneDivSqrtR2pR1=1.;
        OneDivR2pR1    =1.;
        R1DivR1pR2     =0.;
    }
    assert(basis.d==2);
    bool is_realline_basis = flens::IsSame<Basis<T,Primal,R,CDF>, Basis1D>::value;
    assert(is_realline_basis || (R1>0 && R2>0.));
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, BlackScholes, Basis1D>::operator()(XType xtype, int j, int k) const
{
    if (!excessToPayoff) return 0.;
    T sigma    = processparameters.sigma;
    T strike   = optionparameters.strike;
    return 0.5*sigma*sigma* strike * OneDivSqrtR2pR1 * basis.generator(xtype)(R1DivR1pR2,j,k,0);
}


template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, BlackScholes, Basis1D>::operator()(T time, XType xtype, int j, int k) const
{
    return this->operator()(xtype, j, k);
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, BlackScholes, Basis1D>::operator()(const Index1D &lambda) const
{
    return this->operator()(lambda.xtype, lambda.j, lambda.k);
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, BlackScholes, Basis1D>::operator()(T time, const Index1D &lambda) const
{
    return this->operator()(lambda.xtype, lambda.j, lambda.k);
}

}   // namespace lawa
