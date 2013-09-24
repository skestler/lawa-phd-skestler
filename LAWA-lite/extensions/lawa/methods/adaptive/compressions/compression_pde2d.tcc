namespace lawa {

template <typename T, typename Basis2D>
CompressionPDE2D<T,Basis2D>::CompressionPDE2D(const Basis2D &_basis, bool _levelthresh, short _J)
    : basis(_basis), levelthresh(_levelthresh), J(_J),
      s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100), jmax_y(-30)
{
    assert(basis.first.d==basis.second.d);
}

template <typename T, typename Basis2D>
void
CompressionPDE2D<T,Basis2D>::setParameters(const IndexSet<Index2D> &LambdaRow) {
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
CompressionPDE2D<T,Basis2D>::SparsityPattern(const Index2D &lambda_col,
                                             const IndexSet<Index2D> &LambdaRow)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;

    IndexSet<Index2D> LambdaRowSparse;
    IndexSet<Index1D> Lambda_x =
               lambdaTilde1d_PDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x, jmax_x, false);
    IndexSet<Index1D> Lambda_y =
               lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y, jmax_y, false);

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

/*
template <typename T, typename Basis2D>
IndexSet<Index2D>
CompressionPDE2D<T,Basis2D>::SparsityPattern(const Index2D &lambda_col, int jmin_x, int jmin_y,
                                             int s_tilde, int deriv_x, int deriv_y) {
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    IndexSet<Index2D> ret(basis.first.d,basis.first.d_);

    T factor_x = basis.first.d -  2+1.5-deriv_x;
    T factor_y = basis.second.d - 2+1.5-deriv_y;

    int level_bound_x = round((1./factor_x)*s_tilde);
    int level_bound_y = round((1./factor_y)*s_tilde);

    for (int s_tilde_x=0; s_tilde_x <= level_bound_x; ++s_tilde_x) {
        for (int s_tilde_y=0; s_tilde_y <= level_bound_y; ++s_tilde_y) {
            if (factor_x*s_tilde_x + factor_y*s_tilde_y <= s_tilde) {
                IndexSet<Index1D> Lambda_x =
                        lambdaTilde1d_PDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x,
                                          lambda_col.index1.j+s_tilde_x, false);
                IndexSet<Index1D> Lambda_y =
                        lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y,
                                          lambda_col.index2.j+s_tilde_y, false);


                for (set1d_const_it lambda_x=Lambda_x.begin(); lambda_x!=Lambda_x.end(); ++lambda_x) {
                    for (set1d_const_it lambda_y=Lambda_y.begin(); lambda_y!=Lambda_y.end(); ++lambda_y) {
                        if (deriv_x!=0 || deriv_y!=0) {
                            ret.insert(Index2D(*lambda_x,*lambda_y));
                        }
                        else {
                            int tmp =   std::max(lambda_col.index1.j,lambda_col.index2.j)
                                      + std::max((*lambda_x).j,(*lambda_y).j);
                            if (tmp <= s_tilde) {
                                ret.insert(Index2D(*lambda_x,*lambda_y));
                            }
                        }
                    }
                }
            }
        }
    }

    return ret;
}
*/

}    //namespace lawa

