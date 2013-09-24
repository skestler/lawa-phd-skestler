namespace lawa {

template <typename T, typename Basis2D>
IdentityOperator2D<T, Basis2D>::IdentityOperator2D(const Basis2D &_basis)
: basis(_basis),
  integral_x(_basis.first, _basis.first), integral_y(_basis.second, _basis.second)
{
}

template <typename T, typename Basis2D>
T
IdentityOperator2D<T, Basis2D>::operator()(XType row_xtype_x, int j1_x, long k1_x,
                                           XType row_xtype_y, int j1_y, long k1_y,
                                           XType col_xtype_x, int j2_x, long k2_x,
                                           XType col_xtype_y, int j2_y, long k2_y) const
{   
    // (v1 * u1) * (v2 * u2)
    return integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0) 
         * integral_y(j1_y, k1_y, row_xtype_y, 0, j2_y, k2_y, col_xtype_y, 0);
}

template <typename T, typename Basis2D>
T
IdentityOperator2D<T, Basis2D>::operator()(const Index2D &row_index, const Index2D &col_index) const
{
    return HelmholtzOperator2D<T, Basis2D>::operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                                                       row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                                                       col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                                                       col_index.index2.xtype, col_index.index2.j, col_index.index2.k);
}

}   //namespace lawa

