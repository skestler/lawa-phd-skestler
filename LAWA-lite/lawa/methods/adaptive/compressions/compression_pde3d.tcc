namespace lawa {

template <typename T, typename Basis3D>
CompressionPDE3D<T,Basis3D>::CompressionPDE3D(const Basis3D &_basis)
    : basis(_basis), s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100),
      jmax_y(-30), s_tilde_z(-1), jmin_z(100), jmax_z(-30)
{
}

template <typename T, typename Basis3D>
void
CompressionPDE3D<T,Basis3D>::setParameters(const IndexSet<Index3D> &LambdaRow)
{
    typedef typename IndexSet<Index3D>::const_iterator set3d_const_it;
    jmin_x = 100, jmax_x=-30, jmin_y = 100, jmax_y=-30, jmin_z = 100, jmax_z=-30;
    for (set3d_const_it lambda_col=LambdaRow.begin(); lambda_col!=LambdaRow.end(); ++lambda_col) {
        jmin_x = std::min(jmin_x,(*lambda_col).index1.j);
        jmax_x = std::max(jmax_x,(*lambda_col).index1.j);
        jmin_y = std::min(jmin_y,(*lambda_col).index2.j);
        jmax_y = std::max(jmax_y,(*lambda_col).index2.j);
        jmin_z = std::min(jmin_z,(*lambda_col).index3.j);
        jmax_z = std::max(jmax_z,(*lambda_col).index3.j);
    }
    s_tilde_x = jmax_x-jmin_x;
    s_tilde_y = jmax_y-jmin_y;
    s_tilde_z = jmax_z-jmin_z;
}

template <typename T, typename Basis3D>
IndexSet<Index3D>
CompressionPDE3D<T,Basis3D>::SparsityPattern(const Index3D &lambda_col,
                                             const IndexSet<Index3D> &LambdaRow)
{
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
    typedef typename IndexSet<Index3D>::const_iterator set3d_const_it;

    set3d_const_it LambdaRow_end = LambdaRow.end();
    IndexSet<Index3D> LambdaRowSparse;

    IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first, s_tilde_x,
                                                   jmin_x, jmax_x, false);
    IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y,
                                                   jmin_y, jmax_y, false);
    IndexSet<Index1D> Lambda_z = lambdaTilde1d_PDE(lambda_col.index3, basis.third, s_tilde_z,
                                                   jmin_z, jmax_z, false);

    for (set3d_const_it it=LambdaRow.begin(); it!=LambdaRow.end(); ++it) {
        if (Lambda_x.count((*it).index1) >0 ) {
            if (Lambda_y.count((*it).index2) >0 ) {
                if (Lambda_z.count((*it).index3) >0 ) {
                    LambdaRowSparse.insert(Index3D((*it).index1,(*it).index2,(*it).index3));
                }
            }
        }
    }
/*
    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
        for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
            for (set1d_const_it lambda_z = Lambda_z.begin(); lambda_z != Lambda_z.end(); ++lambda_z) {
                Index3D index3d(*lambda_x,*lambda_y,*lambda_z);
                set3d_const_it it = LambdaRow.find(index3d);
                if (it != LambdaRow_end) {
                    LambdaRowSparse.insert(index3d);
                }
                //if (LambdaRow.count(index3d) > 0) {
                //    LambdaRowSparse.insert(index3d);
                //}
            }
        }
    }
*/
    return LambdaRowSparse;
}

}    //namespace lawa

