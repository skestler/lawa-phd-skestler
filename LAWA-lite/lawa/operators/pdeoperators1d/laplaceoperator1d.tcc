namespace lawa {

template <typename T, typename Basis>
LaplaceOperator1D<T,Basis>::LaplaceOperator1D(const Basis& _basis)
    : basis(_basis), integral(_basis, _basis)
{
}

template <typename T, typename Basis>
T
LaplaceOperator1D<T,Basis>::operator()(XType xtype1, int j1, long k1,
                                       XType xtype2, int j2, long k2) const
{   
    // v_x * u_x
    return integral(j1, k1, xtype1, 1, j2, k2, xtype2, 1);
}

template <typename T, typename Basis>
T
LaplaceOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return LaplaceOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                   col_index.xtype, col_index.j, col_index.k);
}

}    //namespace lawa
