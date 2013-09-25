namespace lawa {

template <typename T>
T
StabilizingExpFunction<T>::lambda;

template <typename T>
flens::DenseVector<Array<T> >
StabilizingExpFunction<T>::sing_pts;


template <typename T>
void
StabilizingExpFunction<T>::setlambda(T _lambda)
{
    lambda = _lambda;
}

template <typename T>
T
StabilizingExpFunction<T>::ExpPowMlambda(T t)
{
    return std::exp(-lambda*t);
}


template <typename T, typename Basis>
RHSEuropeanPut<T,Basis>::RHSEuropeanPut(const Basis &_basis, T (*f)(T), DenseVector<Array<T> > &sing_pts,
                                        T _C, T _G, T _M, T _Y, T  _r, T _K, T _maturity, T _sigma)
: basis(_basis), stabilizing_exp_function(f,sing_pts),
  phi_t(basis.first.mra), psi_t(basis.first),
  phi_x(basis.second.mra), d_phi_x(basis.second.mra,1), psi_x(basis.second), d_psi_x(basis.second,1),
  integral_sff(phi_t,stabilizing_exp_function), integral_wf(psi_t,stabilizing_exp_function),
  cgmy(_C,_G,_M,_Y), r(_r), K(_K), maturity(_maturity), sigma(_sigma),
  cgmy_operator(basis.second, 0., 0., 0., _C, _G, _M, _Y)
{
    assert(basis.second.d==2);
    integral_sff.quadrature.setOrder(10);
    integral_wf.quadrature.setOrder(10);
}

template <typename T, typename Basis>
T
RHSEuropeanPut<T,Basis>::operator()(XType xtype, int j, int k) const
{
    T ret = 0.;
    GeMatrix<FullStorage<T,ColMajor> > phi_mu_deltas;

    if (xtype==XBSpline) {
        if (sigma != 0) ret += 0.5*sigma*sigma* K * phi_x(0.,j,k);

        if (cgmy.C != 0) {

            T step = pow2i<T>(-j-5);
            ret += K*(cgmy.ExpXmOne_k1_pos+cgmy.ExpXmOne_k1_neg) * phi_x(0,j,k);
            ret -= K*(cgmy.ExpXmOne_k2_neg * d_phi_x(+step,j,k));
            ret -= K*(cgmy.ExpXmOne_k2_pos * d_phi_x(-step,j,k));

            phi_mu_deltas=cgmy_operator.computeDeltas(xtype,j,k);
        }
    }
    else {
        if (sigma != 0) ret += 0.5*sigma*sigma* K * psi_x(0.,j,k);

        if (cgmy.C != 0) {
            T step = pow2i<T>(-j-5);
            ret += K*(cgmy.ExpXmOne_k1_pos+cgmy.ExpXmOne_k1_neg) * psi_x(0,j,k);
            ret -= K*(cgmy.ExpXmOne_k2_neg * d_psi_x(+step,j,k));
            ret -= K*(cgmy.ExpXmOne_k2_pos * d_psi_x(-step,j,k));

            phi_mu_deltas=cgmy_operator.computeDeltas(xtype,j,k);
        }
    }
    if (cgmy.C!=0) {
        for (int i=phi_mu_deltas.rows().firstIndex(); i<=phi_mu_deltas.rows().lastIndex(); ++i) {
            T x = phi_mu_deltas(i,1);
            T c = phi_mu_deltas(i,2);
            if (x != 0) {
                ret += c*( K*exp(x)*cgmy.ThirdTailFirstExpMomIntegral(-x) - K*cgmy.ThirdTailIntegral(-x));
            }
        }
    }
    return ret;

}

template <typename T, typename Basis>
T
RHSEuropeanPut<T,Basis>::operator()(XType xtype_t, int j_t, int k_t, XType xtype_x, int j_x, int k_x) const
{
    T val_t = 0.;
    if (xtype_t==XBSpline) {
        val_t = integral_sff(j_t, k_t);
    }
    else {
        val_t = integral_wf(j_t, k_t);
    }
    if (fabs(val_t)>0) return val_t * operator()(xtype_x,j_x,k_x);
    else               return 0.;
}

template <typename T, typename Basis>
T
RHSEuropeanPut<T,Basis>::operator()(const Index1D &lambda) const
{
    return operator()(lambda.xtype, lambda.j, lambda.k);
}

template <typename T, typename Basis>
T
RHSEuropeanPut<T,Basis>::operator()(const Index2D &lambda) const
{
    return operator()(lambda.index1.xtype, lambda.index1.j, lambda.index1.k,
                      lambda.index2.xtype, lambda.index2.j, lambda.index2.k);
}

}    //namespace lawa
