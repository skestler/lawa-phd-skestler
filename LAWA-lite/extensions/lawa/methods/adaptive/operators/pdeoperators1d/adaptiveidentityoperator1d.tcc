namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
AdaptiveIdentityOperator1D<T,Side,Domain,Cons>::AdaptiveIdentityOperator1D(const Basis1D &_basis1d)
: basis1d(_basis1d),
  compression1d(basis1d), identity_op1d(basis1d), prec1d(),
  data(identity_op1d, prec1d, compression1d)
{

}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveIdentityOperator1D<T,Side,Domain,Cons>::operator()(const Index1D &row_index,
                                                           const Index1D &col_index)
{
    if (Domain==R) {
        int min_j = std::min((int)row_index.j, (int)col_index.j);
        Index1D tmp_row_index(row_index.j-min_j,row_index.k,row_index.xtype);
        Index1D tmp_col_index(col_index.j-min_j,col_index.k,col_index.xtype);

        return data(tmp_row_index, tmp_col_index);
    }
    else {
        return data(row_index, col_index);
    }
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveIdentityOperator1D<T,Side,Domain,Cons>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    //return identity_op1d(xtype_row,j_row,k_row,xtype_col,j_col,k_col);
    return this->operator()(row_index,col_index);
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
void
AdaptiveIdentityOperator1D<T,Side,Domain,Cons>::toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow,
                                                                    const IndexSet<Index1D>& LambdaCol,
                                                                    SparseMatrixT &A_flens, int /*J*/,
                                                                    bool /*useLinearIndex*/)
{
    this->data.toFlensSparseMatrix(LambdaRow, LambdaCol, A_flens,-1);
}

}   // namespace lawa
