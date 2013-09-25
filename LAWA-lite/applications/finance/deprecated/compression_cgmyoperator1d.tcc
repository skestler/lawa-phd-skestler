namespace lawa {

template <typename T>
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::CompressionCGMYOperator1D(const Basis<T,Primal,R,CDF> &_basis, T Y)
    : basis(_basis), jmin(100), jmax(-30), compr_c(1.), psi(basis)
{
    compr_c = T(2*basis.d)/T(2*basis.d_+Y);
}

template <typename T>
void
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::setParameters(const IndexSet<Index1D> &LambdaRow) {
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    jmin = 100;
    jmax = -30;
    for (set1d_const_it lambda_row = LambdaRow.begin(); lambda_row != LambdaRow.end(); ++lambda_row) {
        jmin = std::min(jmin,(*lambda_row).j);
        jmax = std::max(jmax,(*lambda_row).j);
    }
}

template <typename T>
IndexSet<Index1D>
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow) {
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

    IndexSet<Index1D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
    if (lambda_col.xtype == XBSpline) {
        return LambdaRow;
    }
    for (set1d_const_it lambda = LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
        if ((*lambda).xtype == XBSpline) {
            LambdaRowSparse.insert(*lambda);
        }
        else {
            T delta = 3*compr_c*std::max(pow2i<T>(-std::min(lambda_col.j,(*lambda).j)),
                                         pow2i<T>(-(jmax-1) + compr_alpha*(2*(jmax-1)-lambda_col.j-(*lambda).j)));
            if (lawa::distance(psi.support(lambda_col.j,lambda_col.k),psi.support((*lambda).j,(*lambda).k)) < delta) {
                LambdaRowSparse.insert(*lambda);
            }
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa
