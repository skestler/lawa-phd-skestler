namespace lawa{

template <typename T, typename Basis>
SpaceTimeInitialCondition1D<T, Basis>::SpaceTimeInitialCondition1D(const Basis& _basis)
    : basis(_basis), integral_x(_basis.second, _basis.second)
{
}

template <typename T, typename Basis>
T
SpaceTimeInitialCondition1D<T, Basis>::operator()(XType row_xtype_x, int j1_x, long k1_x,
                                                  XType col_xtype_t, int j2_t, long k2_t,
                                                  XType col_xtype_x, int j2_x, long k2_x) const
{
    T factor = basis.first.generator(col_xtype_t)(0, j2_t, k2_t, 0);
    
    // u1(0) * Integral(v2 * u2)
    return factor * integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0);

}

template <typename T, typename Basis>
T
SpaceTimeInitialCondition1D<T, Basis>::operator()(const Index1D &row_index,
                                                  const Index2D &col_index) const
{
    return operator()(row_index.xtype, row_index.j, row_index.k,
                      col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                      col_index.index2.xtype, col_index.index2.j, col_index.index2.k);

}

} // namespace lawa
