namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
AdaptiveLaplaceOperator1D<T,Side,Domain,Cons>::AdaptiveLaplaceOperator1D(const Basis1D &_basis1d)
: basis1d(_basis1d),
  compression1d(basis1d), laplace_op1d(basis1d), prec1d(),
  data(laplace_op1d, prec1d, compression1d)
{

}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveLaplaceOperator1D<T,Side,Domain,Cons>::operator()(const Index1D &row_index,
                                                          const Index1D &col_index)
{
    if (Domain==R) {
        int min_j = std::min((int)row_index.j, (int)col_index.j);
        T scaling_factor = pow2i<T>(2*min_j);
        Index1D tmp_row_index(row_index.j-min_j,row_index.k,row_index.xtype);
        Index1D tmp_col_index(col_index.j-min_j,col_index.k,col_index.xtype);
        return scaling_factor * data(tmp_row_index, tmp_col_index);
    }
    else {
        return data(row_index, col_index);
    }
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveLaplaceOperator1D<T,Side,Domain,Cons>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    //return laplace_op1d(xtype_row,j_row,k_row,xtype_col,j_col,k_col);
    return this->operator()(row_index,col_index);
}

}   // namespace lawa
