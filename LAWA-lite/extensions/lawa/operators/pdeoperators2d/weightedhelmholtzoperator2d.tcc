namespace lawa {

template <typename T, typename Basis>
WeightedHelmholtzOperator2D<T, Basis>::WeightedHelmholtzOperator2D(const Basis& _basis, const T _c,
                                                               Function<T> weightFct_x, 
                                                               Function<T> weightFct_y)
    : basis(_basis), c(_c), W_x(weightFct_x), W_y(weightFct_y),
      integral_x(W_x, basis.first, basis.first), integral_y(W_y, basis.second, basis.second)
{}


template <typename T, typename Basis>
T
WeightedHelmholtzOperator2D<T, Basis>::operator()(XType row_xtype_x, int j1_x, long k1_x,
                                                  XType row_xtype_y, int j1_y, long k1_y,
                                                  XType col_xtype_x, int j2_x, long k2_x,
                                                  XType col_xtype_y, int j2_y, long k2_y) const
{
    return 	integral_x(j1_x, k1_x, row_xtype_x, 1, j2_x, k2_x, col_xtype_x, 1) 
          * integral_y(j1_y, k1_y, row_xtype_y, 0, j2_y, k2_y, col_xtype_y, 0)
         +  integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0) 
          * integral_y(j1_y, k1_y, row_xtype_y, 1, j2_y, k2_y, col_xtype_y, 1)
      + c * integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0)
          * integral_y(j1_y, k1_y, row_xtype_y, 0, j2_y, k2_y, col_xtype_y, 0);
}

template <typename T, typename Basis>
T
WeightedHelmholtzOperator2D<T, Basis>::operator()(const Index2D &row_index, const Index2D &col_index) const
{
      return this->operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                              row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                              col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                              col_index.index2.xtype, col_index.index2.j, col_index.index2.k);
}

template <typename T, typename Basis>
T
WeightedHelmholtzOperator2D<T, Basis>::operator()(const Index2D &row_index, const Index2D &col_index)
{
      return this->operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                              row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                              col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                              col_index.index2.xtype, col_index.index2.j, col_index.index2.k);
}

} // namespace lawa
