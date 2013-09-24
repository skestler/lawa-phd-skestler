namespace lawa {

template <typename T, typename Basis, typename Preconditioner>
WeightedAdaptiveHelmholtzOperator2D<T, Basis, Preconditioner>::WeightedAdaptiveHelmholtzOperator2D
           (const Basis &_basis, const T _c, Function<T> weightFct_x, Function<T> weightFct_y,
            const Preconditioner &_Prec)
: basis(_basis), c(_c), Prec(_Prec),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  op_identity_x(basis.first, weightFct_x), op_identity_y(basis.second, weightFct_y),
  op_laplace_x(basis.first, weightFct_x), op_laplace_y(basis.second, weightFct_y),
  Prec1D(),
  data_identity_x(op_identity_x,   Prec1D, compression_1d_x),
  data_identity_y(op_identity_y,   Prec1D, compression_1d_y),
  data_laplace_x(op_laplace_x,     Prec1D, compression_1d_x),
  data_laplace_y(op_laplace_y,     Prec1D, compression_1d_y),
  P_data()
{

}

template <typename T, typename Basis, typename Preconditioner>
T
WeightedAdaptiveHelmholtzOperator2D<T, Basis, Preconditioner>::operator()(const Index2D &row_index,
                                                                    const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    if (!flens::IsSame<NoPreconditioner2D, Preconditioner>::value) {
        // Left precondioning:
        const_coeff_it it_row_index   = P_data.find(row_index);
        //  Entry has already been computed:
        if (it_row_index != P_data.end()) {
            prec *= (*it_row_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = Prec(row_index);
            P_data[row_index] = tmp;
            prec *= tmp;
        }
        // Right precondioning:
        const_coeff_it it_col_index   = P_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = Prec(col_index);
            P_data[col_index] = tmp;
            prec *= tmp;
        }
    }

    T dd_x = data_laplace_x(row_index.index1,col_index.index1);
    T id_x = data_identity_x(row_index.index1,col_index.index1);
    T dd_y = data_laplace_y(row_index.index2,col_index.index2);
    T id_y = data_identity_y(row_index.index2,col_index.index2);

    return prec*(dd_x*id_y + id_x*dd_y + c*id_x*id_y  );
}

template <typename T, typename Basis, typename Preconditioner>
T
WeightedAdaptiveHelmholtzOperator2D<T, Basis, Preconditioner>::prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    if (!flens::IsSame<NoPreconditioner2D, Preconditioner>::value) {
        const_coeff_it it_index   = P_data.find(index);
        if (it_index != P_data.end()) {
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

template <typename T, typename Basis, typename Preconditioner>
void
WeightedAdaptiveHelmholtzOperator2D<T, Basis, Preconditioner>::clear()
{
    data_laplace_x.clear();
    data_identity_x.clear();
    data_laplace_y.clear();
    data_identity_y.clear();
}


}   // namespace lawa

