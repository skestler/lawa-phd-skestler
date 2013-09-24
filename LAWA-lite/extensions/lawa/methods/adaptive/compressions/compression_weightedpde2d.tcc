namespace lawa {

template <typename T, typename Basis2D>
CompressionWeightedPDE2D<T,Basis2D>::CompressionWeightedPDE2D(const Basis2D &_basis, bool _levelthresh, short _J)
    : basis(_basis), levelthresh(_levelthresh), J(_J),
      s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100), jmax_y(-30)
{
    assert(basis.first.d==basis.second.d);
}

template <typename T, typename Basis2D>
void
CompressionWeightedPDE2D<T,Basis2D>::setParameters(const IndexSet<Index2D> &LambdaRow) {
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;
    jmin_x = 100, jmax_x=-30, jmin_y = 100, jmax_y=-30;
    for (set2d_const_it lambda_col=LambdaRow.begin(); lambda_col!=LambdaRow.end(); ++lambda_col) {
        jmin_x = std::min(jmin_x,(*lambda_col).index1.j);
        jmax_x = std::max(jmax_x,(*lambda_col).index1.j);
        jmin_y = std::min(jmin_y,(*lambda_col).index2.j);
        jmax_y = std::max(jmax_y,(*lambda_col).index2.j);
    }
    s_tilde_x = jmax_x-jmin_x;
    s_tilde_y = jmax_y-jmin_y;
}

template <typename T, typename Basis2D>
IndexSet<Index2D>
CompressionWeightedPDE2D<T,Basis2D>::SparsityPattern(const Index2D &lambda_col,
                                             const IndexSet<Index2D> &LambdaRow)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;

    IndexSet<Index2D> LambdaRowSparse;
    IndexSet<Index1D> Lambda_x =
               lambdaTilde1d_WeightedPDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x, jmax_x);
    IndexSet<Index1D> Lambda_y =
               lambdaTilde1d_WeightedPDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y, jmax_y);

    //int level_thresh_bound = std::min(J,std::max(s_tilde_x,s_tilde_y));

    for (set2d_const_it lambda=LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
        short level_diff =   fabs((*lambda).index1.j-lambda_col.index1.j)
                           + fabs((*lambda).index2.j-lambda_col.index2.j);
        if ( (levelthresh) && ((0.5+basis.first.d-2)*level_diff > J) ) {
            continue;
        }
        if ((Lambda_x.count((*lambda).index1)>0) && (Lambda_y.count((*lambda).index2)>0))  {
            LambdaRowSparse.insert(*lambda);
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa

