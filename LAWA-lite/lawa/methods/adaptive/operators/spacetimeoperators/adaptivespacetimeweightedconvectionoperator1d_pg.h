/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
 
#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_SPACETIMEOPERATORS_ADAPTIVESPACETIMEWEIGHTEDCONVECTIONOPERATOR1D_PG_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_SPACETIMEOPERATORS_ADAPTIVESPACETIMEWEIGHTEDCONVECTIONOPERATOR1D_PG_H 1
 
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/compressions/compression_pde1d.h>
#include <lawa/methods/adaptive/compressions/compression_pde2d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptiveoperator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedconvectionoperator1d_pg.h>
#include <lawa/operators/spacetimeoperators/spacetimeoperators.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/preconditioners/spacetimepreconditioners/spacetimepreconditioners.h>
 
namespace lawa {

/* Space-Time Weighted Convection Operator
 *
 *  a(v,u) =  convection * Integral(v1 * u1)    * Integral( w * v2 * u2_x)
 *
 *  Template Parameters:
 *      LeftPrec2D :        left preconditioner
 *      RightPrec2D:        right preconditioner
 *      InitialCondition:   operator for initial condition, can be NoInitialCondition
 */
template <typename T, typename TrialBasis, typename TestBasis,
          typename TrialPrec, typename TestPrec, typename InitialCondition>
struct AdaptiveSpaceTimeWeightedConvectionOperator1D_PG : public AdaptiveOperator2D<T> {
    
    typedef flens::SparseGeMatrix<CRS<T,CRS_General> >                  SparseMatrixT;

    typedef typename TrialBasis::FirstBasisType    TrialBasis_t;
    typedef typename TrialBasis::SecondBasisType   TrialBasis_x;

    typedef typename TestBasis::FirstBasisType    TestBasis_t;
    typedef typename TestBasis::SecondBasisType   TestBasis_x;
    
    typedef NoCompression<T, Index1D, TrialBasis_t>   Compression1D_t;
    typedef NoCompression<T, Index1D, TrialBasis_x>   Compression1D_x;
    typedef NoCompression<T, Index2D, TrialBasis>     Compression2D;
    
    typedef NoPreconditioner<T,Index1D>         NoPreconditioner1D;
    
    typedef IdentityOperator1D_PG<T, TrialBasis_t, TestBasis_t>      				IdentityOperator_t;
    typedef WeightedConvectionOperator1D_PG<T, TrialBasis_x, TestBasis_x, Gauss>    ConvectionOperator_x;

    typedef MapMatrix<T, Index1D, IdentityOperator_t, 
                               Compression1D_t, NoPreconditioner1D>   DataIdentity_t;  
    typedef MapMatrix<T, Index1D, ConvectionOperator_x, 
                               Compression1D_x, NoPreconditioner1D>   DataConvection_x;
                               
    AdaptiveSpaceTimeWeightedConvectionOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
            TrialPrec& _trialprec, TestPrec& _testprec, Function<T> weightFct, T _convection = 1);
    
    AdaptiveSpaceTimeWeightedConvectionOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
            TrialPrec& _trialprec, TestPrec& _testprec, InitialCondition& _init_cond, Function<T> weightFct, T _convection = 1);
                                    
    // call of p_left * a_operator * p_right
    T
    operator()(const Index2D &row_index, const Index2D &col_index);
    
    //call of op_initcond * p_right
    T
    operator()(const Index1D &row_index, const Index2D &col_index);

    Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow,
                        const IndexSet<Index2D> &LambdaCol, SparseMatrixT &A, T tol, bool useLinearIndex=false);

    void
    clear();
    
    const TrialBasis&   trialbasis;
    const TestBasis&    testbasis;
    const T             convection;
    
    Compression1D_t     compression_1d_t;
    Compression1D_x     compression_1d_x;
    Compression2D       compression;
    
    Coefficients<Lexicographical,T,Index2D> P_trial_data;
    Coefficients<Lexicographical,T,Index2D> P_test_data;
    
    const TrialPrec&   trialprec;
    const TestPrec&    testprec;
    NoPreconditioner1D  noprec;
    
    const IdentityOperator_t    op_identity_t;
    const ConvectionOperator_x  op_convection_x;
    
    const NoInitialCondition    op_noinitcond;
    const InitialCondition&     op_initcond;
    
    DataIdentity_t      data_identity_t;
    DataConvection_x    data_convection_x;
    
};    
      
} // namespace lawa

#include <lawa/methods/adaptive/operators/spacetimeoperators/adaptivespacetimeweightedconvectionoperator1d_pg.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_SPACETIMEOPERATORS_ADAPTIVESPACETIMEWEIGHTEDCONVECTIONOPERATOR1D_PG_H
 
 

