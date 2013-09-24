namespace lawa {

template <typename T, typename Basis2D, typename Preconditioner>
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::AdaptiveHelmholtzOperator2D
                                                         (const Basis2D &_basis, T _c,
                                                          const Preconditioner &_Prec,
                                                          T /*thresh*/, int /*NumOfCols*/,
                                                          int /*NumOfRows*/)
: basis(_basis), c(_c), Prec(_Prec),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  op_identity_x(basis.first), op_identity_y(basis.second),
  op_laplace_x(basis.first), op_laplace_y(basis.second),
  Prec1D(),
  data_identity_x(op_identity_x, Prec1D, compression_1d_x),
  data_identity_y(op_identity_y, Prec1D, compression_1d_y),
  data_laplace_x(op_laplace_x,   Prec1D, compression_1d_x),
  data_laplace_y(op_laplace_y,   Prec1D, compression_1d_y),
  P_data()
{

}

template <typename T, typename Basis2D, typename Preconditioner>
T
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::operator()(const Index2D &row_index,
                                                                    const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;

    const_coeff_it it_P_end       = P_data.end();
    const_coeff_it it_row_index   = P_data.find(row_index);
    if (it_row_index != it_P_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp;
        if (flens::IsSame<DiagonalHelmholtzPreconditioner2D, Preconditioner>::value) {
            T prec_dd_x = data_laplace_x(row_index.index1,row_index.index1);
            T prec_id_x = data_identity_x(row_index.index1,row_index.index1);
            T prec_id_y = data_identity_y(row_index.index2,row_index.index2);
            T prec_dd_y = data_laplace_y(row_index.index2,row_index.index2);
            tmp = 1./std::sqrt(fabs(prec_dd_x*prec_id_y + prec_id_x*prec_dd_y
                                    + c*prec_id_x*prec_id_y ));
        }
        else {
            tmp = Prec(row_index);
        }
        P_data[row_index] = tmp;
        prec *= tmp;
    }
    it_P_end       = P_data.end();
    const_coeff_it it_col_index   = P_data.find(col_index);
    if (it_col_index != it_P_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp;
        if (flens::IsSame<DiagonalHelmholtzPreconditioner2D, Preconditioner>::value) {
            T prec_dd_x = data_laplace_x(col_index.index1,col_index.index1);
            T prec_id_x = data_identity_x(col_index.index1,col_index.index1);
            T prec_id_y = data_identity_y(col_index.index2,col_index.index2);
            T prec_dd_y = data_laplace_y(col_index.index2,col_index.index2);
            tmp = 1./std::sqrt(fabs(prec_dd_x*prec_id_y + prec_id_x*prec_dd_y
                                    + c*prec_id_x*prec_id_y ));
        }
        else {
            tmp = Prec(col_index);
        }
        P_data[col_index] = tmp;
        prec *= tmp;
    }

    T dd_x = data_laplace_x(row_index.index1,col_index.index1);
    T id_x = data_identity_x(row_index.index1,col_index.index1);

    T id_y, dd_y;
    if ( (flens::IsSame<typename Basis2D::FirstBasisType, Basis<T,Primal,R,CDF> >::value) &&
         (flens::IsSame<typename Basis2D::SecondBasisType, Basis<T,Primal,R,CDF> >::value)    ) {
        if (row_index.index2.xtype==XWavelet && col_index.index2.xtype==XWavelet) {
            id_y = data_identity_x(row_index.index2,col_index.index2);
            dd_y = data_laplace_x(row_index.index2,col_index.index2);
        }
        else {
            id_y = data_identity_y(row_index.index2,col_index.index2);
            dd_y = data_laplace_y(row_index.index2,col_index.index2);
        }
    }
    else {
        dd_y = data_laplace_y(row_index.index2,col_index.index2);
        id_y = data_identity_y(row_index.index2,col_index.index2);
    }

    return prec*(dd_x*id_y + id_x*dd_y + c*id_x*id_y  );
}

template <typename T, typename Basis2D, typename Preconditioner>
T
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    if (!flens::IsSame<NoPreconditioner2D, Preconditioner>::value) {
        const_coeff_it it_P_end       = P_data.end();
        const_coeff_it it_index   = P_data.find(index);
        if (it_index != it_P_end) {
            prec *= (*it_index).second;
        }
        else {
            T tmp = Prec(index);
            P_data[index] = tmp;
            prec *= tmp;
        }
    }
    return prec;
}

template <typename T, typename Basis2D, typename Preconditioner>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::mv(const IndexSet<Index2D> &LambdaRow,
                                                const Coefficients<Lexicographical,T,Index2D> &v)
{
    return lawa::mv_sparse(LambdaRow, (*this), v);
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::toFlensSparseMatrix
                                                         (const IndexSet<Index2D>& LambdaRow,
                                                          const IndexSet<Index2D>& LambdaCol,
                                                          SparseMatrixT &A_flens, int /*J*/,
                                                          bool useLinearIndex)
{
    if (!useLinearIndex) {
        lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
    }
    else {
        std::cerr << "AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::toFlensSparseMatrix "
                  << "not implemented for useLinearIndex=true." << std::endl;
        assert(0);
        exit(1);
    }
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::toFlensSparseMatrix(
                                                         const IndexSet<Index2D>& LambdaRow,
                                                         const IndexSet<Index2D>& LambdaCol,
                                                         SparseMatrixT &A_flens, T /*eps*/,
                                                         bool /*useLinearIndex*/)
{
    //toFlensSparseMatrix<T,Index2D,AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner> >(*this,LambdaRow,LambdaCol,A_flens);
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename Basis2D, typename Preconditioner>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::apply
                                                  (const Coefficients<Lexicographical,T,Index2D> &/*v*/,
                                                   int /*k*/, int /*J*/,
                                                   cxxblas::Transpose /*trans*/)
{
    std::cerr << "Apply not yet implemented for this operator." << std::endl;
    exit(1);
    Coefficients<Lexicographical,T,Index2D> ret;
    return ret;
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::apply
                                                  (const Coefficients<Lexicographical,T,Index2D> &/*v*/,
                                                   T /*eps*/, Coefficients<Lexicographical,T,Index2D> &/*ret*/,
                                                   cxxblas::Transpose /*trans*/)
{
    std::cerr << "Apply not yet implemented for this operator." << std::endl;
    exit(1);
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::apply
                                                  (const Coefficients<Lexicographical,T,Index2D> &/*v*/,
                                                   T /*eps*/, const IndexSet<Index2D> &/*Lambda*/,
                                                   Coefficients<Lexicographical,T,Index2D> &/*ret*/,
                                                   cxxblas::Transpose /*trans*/)
{
    std::cerr << "Apply not yet implemented for this operator." << std::endl;
    exit(1);
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::clear()
{
    data_laplace_x.clear();
    data_identity_x.clear();
    data_laplace_y.clear();
    data_identity_y.clear();
}

}   // namespace lawa

