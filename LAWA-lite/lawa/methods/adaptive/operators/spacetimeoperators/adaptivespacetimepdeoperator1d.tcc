#include <iostream>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>

namespace lawa {

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
    AdaptiveSpaceTimePDEOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                   T _diffusion, T _convection, T _reaction, T _timederivfactor)
    : basis(_basis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      timederivfactor(_timederivfactor),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression(_basis),
      P_left_data(), P_right_data(), p_left(_p_left), p_right(_p_right), noprec(),
      op_identity_t(_basis.first), op_identity_x(_basis.second), op_convection_t(_basis.first),
      op_convection_x(_basis.second), op_laplace_x(_basis.second), 
      op_noinitcond(), op_initcond(op_noinitcond),
      data_identity_t(op_identity_t,     noprec, compression_1d_t),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t),
      data_convection_x(op_convection_x, noprec, compression_1d_x),
      data_laplace_x(op_laplace_x,       noprec, compression_1d_x)
{
}
    
template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
    AdaptiveSpaceTimePDEOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                   InitialCondition& _init_cond,
                                   T _diffusion, T _convection, T _reaction, T _timederivfactor)
    : basis(_basis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      timederivfactor(_timederivfactor),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression(_basis),
      P_left_data(), P_right_data(), p_left(_p_left), p_right(_p_right), noprec(),
      op_identity_t(_basis.first), op_identity_x(_basis.second), op_convection_t(_basis.first),
      op_convection_x(_basis.second), op_laplace_x(_basis.second), op_noinitcond(), op_initcond(_init_cond),
      data_identity_t(op_identity_t,     noprec, compression_1d_t),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t),
      data_convection_x(op_convection_x, noprec, compression_1d_x),
      data_laplace_x(op_laplace_x,       noprec, compression_1d_x)
{
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
T
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
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
    
    // Calculate convection and reaction term only if constant not 0
    T convection_term = 0.;
    if (convection != 0) {
        convection_term = data_identity_t(row_index.index1,col_index.index1) * data_convection_x(row_index.index2,col_index.index2);
    }
    T reaction_term = 0.;
    if (reaction != 0) {
        reaction_term = data_identity_t(row_index.index1,col_index.index1) * data_identity_x(row_index.index2,col_index.index2);
    }
    
    return prec * ( timederivfactor * data_convection_t(row_index.index1,col_index.index1) * data_identity_x(row_index.index2,col_index.index2) 
                    + diffusion * data_identity_t(row_index.index1,col_index.index1) * data_laplace_x(row_index.index2,col_index.index2)
                    + convection * convection_term 
                    + reaction * reaction_term);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
T
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
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
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
  return lawa::mv(LambdaRow, *this, x);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
void
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol, bool useLinearIndex)
{
    std::cerr << "AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::"
              << "toFlensSparseMatrix not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
void
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
      const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
      cxxblas::Transpose trans)
{
  std::cerr << "AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::"
			<< "apply not implemented." << std::endl;
  assert(0);
  exit(1);
}

template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
void
AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D, RightPrec2D, InitialCondition>::
clear()
{
    data_identity_t.clear();
    data_identity_x.clear();
    data_convection_t.clear();
    data_convection_x.clear();
    data_laplace_x.clear();
}


} // namespace lawa
