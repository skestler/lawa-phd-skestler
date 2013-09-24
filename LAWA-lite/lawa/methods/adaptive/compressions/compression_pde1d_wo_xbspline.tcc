namespace lawa {

template <typename T, typename Basis>
CompressionPDE1D_WO_XBSpline<T,Basis>::CompressionPDE1D_WO_XBSpline(const Basis &_basis)
    : basis(_basis), s_tilde(0), jmin(100), jmax(-30)
{
}

template <typename T, typename Basis>
void
CompressionPDE1D_WO_XBSpline<T,Basis>::setParameters(const IndexSet<Index1D> &Lambda)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    s_tilde = -1;
    jmin = 100;
    jmax = -30;
    for (set1d_const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        jmin = std::min(jmin,(*lambda).j);
        jmax = std::max(jmax,(*lambda).j);
    }
    s_tilde = jmax-jmin;
}

template <typename T, typename Basis>
IndexSet<Index1D>
CompressionPDE1D_WO_XBSpline<T,Basis>::SparsityPattern(const Index1D &lambda,
                                                       const IndexSet<Index1D> &Lambda, int J)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

    IndexSet<Index1D> LambdaRowSparse, Lambda_x;
    if (J==-1) {
        Lambda_x = lambdaTilde1d_PDE_WO_XBSpline(lambda, basis, s_tilde, jmin, jmax);
    }
    else {
        Lambda_x = lambdaTilde1d_PDE_WO_XBSpline(lambda, basis, std::min(int(s_tilde),J), jmin, jmax);
    }

    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
        if (Lambda.count(*lambda_x)>0) {
            LambdaRowSparse.insert(*lambda_x);
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa

