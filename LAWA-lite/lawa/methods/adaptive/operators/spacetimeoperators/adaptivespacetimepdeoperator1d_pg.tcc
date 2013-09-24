#include <iostream>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
AdaptiveSpaceTimePDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
                                   TrialPrec& _trialprec, TestPrec& _testprec,
                                   T _diffusion, T _convection, T _reaction, T _timederivfactor)
    : trialbasis(_trialbasis), testbasis(_testbasis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      timederivfactor(_timederivfactor),
      compression_1d_t(_trialbasis.first), compression_1d_x(_trialbasis.second), compression(_trialbasis),
      P_trial_data(), P_test_data(), trialprec(_trialprec), testprec(_testprec), noprec(),
      op_identity_t(_trialbasis.first, _testbasis.first), op_identity_x(_trialbasis.second, _testbasis.second), 
      op_convection_t(_trialbasis.first, _testbasis.first), op_convection_x(_trialbasis.second, _testbasis.second), 
      op_laplace_x(_trialbasis.second, _testbasis.second), op_noinitcond(), op_initcond(op_noinitcond),
      data_identity_t(op_identity_t,     noprec, compression_1d_t),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t),
      data_convection_x(op_convection_x, noprec, compression_1d_x),
      data_laplace_x(op_laplace_x,       noprec, compression_1d_x)
{
}
    
template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
AdaptiveSpaceTimePDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
                                   TrialPrec& _trialprec, TestPrec& _testprec,
                                   InitialCondition& _init_cond,
                                   T _diffusion, T _convection, T _reaction, T _timederivfactor)
    : trialbasis(_trialbasis), testbasis(_testbasis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      timederivfactor(_timederivfactor),
      compression_1d_t(_trialbasis.first), compression_1d_x(_trialbasis.second), compression(_trialbasis),
      P_trial_data(), P_test_data(), trialprec(_trialprec), testprec(_testprec), noprec(),
      op_identity_t(_trialbasis.first, _testbasis.first), op_identity_x(_trialbasis.second, _testbasis.second), 
      op_convection_t(_trialbasis.first, _testbasis.first), op_convection_x(_trialbasis.second, _testbasis.second), 
      op_laplace_x(_trialbasis.second, _testbasis.second), op_noinitcond(), op_initcond(op_noinitcond),
      data_identity_t(op_identity_t,     noprec, compression_1d_t),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_convection_t(op_convection_t, noprec, compression_1d_t),
      data_convection_x(op_convection_x, noprec, compression_1d_x),
      data_laplace_x(op_laplace_x,       noprec, compression_1d_x)
{
}

template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
T
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TrialPrec>::value) {
        // Left precondioning:
        const_coeff_it it_row_index   = P_test_data.find(row_index);
        //  Entry has already been computed:
        if (it_row_index != P_test_data.end()) {
            prec *= (*it_row_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = testprec(row_index);
            P_test_data[row_index] = tmp;
            prec *= tmp;
        }
    }

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TestPrec>::value) {
        // Right precondioning:
        const_coeff_it it_col_index   = P_trial_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_trial_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = trialprec(col_index);
            P_trial_data[col_index] = tmp;
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

template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
T
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
operator()(const Index1D &row_index, const Index2D &col_index)
{
    if(flens::IsSame<NoInitialCondition, InitialCondition>::value){
        std::cerr << " Operator cannot be called without Initial Condition " << std::endl;
        exit(1);
    }
    
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TestPrec>::value) {
        // Right precondioning:
        const_coeff_it it_col_index   = P_trial_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_trial_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = trialprec(col_index);
            P_trial_data[col_index] = tmp;
            prec *= tmp;
        }
    }

    return prec * op_initcond(row_index,col_index);
}

template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
Coefficients<Lexicographical,T,Index2D>
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
  return lawa::mv(LambdaRow, *this, x);
}


template <typename T, typename TrialBasis, typename TestBasis,
          typename TrialPrec, typename TestPrec, typename InitialCondition>
void
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol, bool useLinearIndex)
{
    std::cerr << "AdaptiveSpaceTimeConvectionOperator1D<T, Basis2D, TrialPrec, TestPrec, InitialCondition>::"
              << "toFlensSparseMatrix not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename TrialBasis, typename TestBasis, 
          typename TrialPrec, typename TestPrec, typename InitialCondition>
void
AdaptiveSpaceTimePDEOperator1D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, InitialCondition>::
clear()
{
    data_identity_t.clear();
    data_identity_x.clear();
    data_convection_t.clear();
    data_convection_x.clear();
    data_laplace_x.clear();
}
    
} // namespace lawa
