namespace lawa {

template <typename T, typename Basis2D>
AdaptiveIdentityOperator2D<T, Basis2D>::AdaptiveIdentityOperator2D(const Basis2D &_basis,
                                                                    T _entrybound,
                                                                    int _NumOfRows, int _NumOfCols)
: basis(_basis), compression_1d_x(basis.first), compression_1d_y(basis.second),
  compression(basis),
  Prec(), op_identity_x(basis.first), op_identity_y(basis.second),
  entrybound(_entrybound), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
  data_identity_x(op_identity_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
  data_identity_y(op_identity_y, compression_1d_y, entrybound, NumOfRows, NumOfCols)
{

}

template <typename T, typename Basis2D>
T
AdaptiveIdentityOperator2D<T, Basis2D>::operator()(const Index2D &row_index,
                                                   const Index2D &col_index)
{
    return  data_identity_x(row_index.index1,col_index.index1)
           *data_identity_y(row_index.index2,col_index.index2);
}

template <typename T, typename Basis2D>
T
AdaptiveIdentityOperator2D<T, Basis2D>::prec(const Index2D &/*index*/)
{
    return 1.;
}

}   //namespace lawa

