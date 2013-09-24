#include <iostream>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
    AdaptiveSpaceTimeTimeDerivOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                   T _timederivfactor)
    : basis(_basis), timederivfactor(_timederivfactor),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression(_basis),
      P_left_data(), P_right_data(), p_left(_p_left), p_right(_p_right), noprec(),
      op_identity_x(_basis.second), op_convection_t(_basis.first),
      op_noinitcond(), op_initcond(op_noinitcond),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t)
{
}
    
template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
    AdaptiveSpaceTimeTimeDerivOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                   InitialCondition& _init_cond,
                                   T _timederivfactor)
    : basis(_basis), timederivfactor(_timederivfactor),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression(_basis),
      P_left_data(), P_right_data(), p_left(_p_left), p_right(_p_right), noprec(),
      op_identity_x(_basis.second), op_convection_t(_basis.first),
      op_noinitcond(), op_initcond(_init_cond),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t)
{
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
T
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, LeftPrec2D>::value) {
        // Left precondioning:
        const_coeff_it it_row_index   = P_left_data.find(row_index);
        //  Entry has already been computed:
        if (it_row_index != P_left_data.end()) {
            prec *= (*it_row_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = p_left(row_index);
            P_left_data[row_index] = tmp;
            prec *= tmp;
        }
    }

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, RightPrec2D>::value) {
        // Right precondioning:
        const_coeff_it it_col_index   = P_right_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_right_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = p_right(col_index);
            P_right_data[col_index] = tmp;
            prec *= tmp;
        }
    }
    
    return prec * timederivfactor * data_convection_t(row_index.index1,col_index.index1) * data_identity_x(row_index.index2,col_index.index2);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
T
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
operator()(const Index1D &row_index, const Index2D &col_index)
{
    if(flens::IsSame<NoInitialCondition, InitialCondition>::value){
        std::cerr << " Operator cannot be called without Initial Condition " << std::endl;
        exit(1);
    }
    
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, RightPrec2D>::value) {
        // Right precondioning:
        const_coeff_it it_col_index   = P_right_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_right_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = p_right(col_index);
            P_right_data[col_index] = tmp;
            prec *= tmp;
        }
    }

    return prec * op_initcond(row_index,col_index);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
Coefficients<Lexicographical,T,Index2D>
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
    std::cerr << "AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::"
              << "mv not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
void
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol, bool useLinearIndex)
{
    std::cerr << "AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::"
              << "toFlensSparseMatrix not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
void
AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
clear()
{
    data_identity_x.clear();
    data_convection_t.clear();
}


} // namespace lawa
