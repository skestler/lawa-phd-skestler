namespace lawa {

template <typename T, typename Basis>
CompressionPDE1D<T,Basis>::CompressionPDE1D(const Basis &_basis)
    : basis(_basis), s_tilde(-1), jmin(100), jmax(-30)
{
}

template <typename T, typename Basis>
void
CompressionPDE1D<T,Basis>::setParameters(const IndexSet<Index1D> &LambdaRow)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

    s_tilde = -1;
    jmin = 100;
    jmax = -30;
    for (set1d_const_it lambda_row=LambdaRow.begin(); lambda_row!=LambdaRow.end(); ++lambda_row) {
        jmin = std::min(jmin,(*lambda_row).j);
        jmax = std::max(jmax,(*lambda_row).j);
    }
    s_tilde = jmax-jmin;
}

template <typename T, typename Basis>
IndexSet<Index1D>
CompressionPDE1D<T,Basis>::SparsityPattern(const Index1D &lambda_col,
                                           const IndexSet<Index1D> &LambdaRow, int J)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

    IndexSet<Index1D> LambdaRowSparse;
    int s = std::max(abs(lambda_col.j-jmin),abs(lambda_col.j-jmax));
    s = std::max(s,int(s_tilde));
    if (J!=-1) s = std::min(s,J);
    //Compression level J>s_tilde as indices corresponding to level differences
    //larger than s_tilde cannot appear in LambdaRow

    IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col, basis, s, jmin, jmax, false);
    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
        if (LambdaRow.count(*lambda_x)>0) {
            LambdaRowSparse.insert(*lambda_x);
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa

