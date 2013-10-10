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
 
#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVESPACETIMEREACTIONOPERATOR1D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVESPACETIMEREACTIONOPERATOR1D_H 1
 
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/compressions/compression_pde1d.h>
#include <lawa/methods/adaptive/compressions/compression_pde2d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptiveoperator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d.h>
#include <lawa/operators/spacetimeoperators/spacetimeoperators.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/preconditioners/spacetimepreconditioners/spacetimepreconditioners.h>
 
namespace lawa {

/* Space-Time Reaction Operator
 *
 *  a(v,u) =  reaction   * Integral(v1 * u1)    * Integral(v2 * u2)
 *
 *  Template Parameters:
 *      LeftPrec2D :        left preconditioner
 *      RightPrec2D:        right preconditioner
 *      InitialCondition:   operator for initial condition, can be NoInitialCondition
 */
template <typename T, typename Basis2D, typename LeftPrec2D, typename RightPrec2D, typename InitialCondition>
struct AdaptiveSpaceTimeReactionOperator1D : public AdaptiveOperator2D<T> {
    
    typedef flens::SparseGeMatrix<CRS<T,CRS_General> >                  SparseMatrixT;

    typedef typename Basis2D::FirstBasisType    Basis_t;
    typedef typename Basis2D::SecondBasisType   Basis_x;
    
    typedef CompressionPDE1D<T, Basis_t>        Compression1D_t;
    typedef CompressionPDE1D<T, Basis_x>        Compression1D_x;  
    typedef CompressionPDE2D<T, Basis2D>        Compression2D;
    
    typedef NoPreconditioner<T,Index1D>         NoPreconditioner1D;
    
    typedef IdentityOperator1D<T, Basis_t>      IdentityOperator_t;
    typedef IdentityOperator1D<T, Basis_x>      IdentityOperator_x;

    typedef MapMatrix<T, Index1D, IdentityOperator_t, 
                               Compression1D_t, NoPreconditioner1D>   DataIdentity_t;
    typedef MapMatrix<T, Index1D, IdentityOperator_x, 
                               Compression1D_x, NoPreconditioner1D>   DataIdentity_x;    
                               
    AdaptiveSpaceTimeReactionOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                    T _reaction = 1.);
    
    AdaptiveSpaceTimeReactionOperator1D(const Basis2D& _basis, LeftPrec2D& _p_left, RightPrec2D& _p_right,
                                    InitialCondition& _init_cond,
                                    T _reaction = 1.);
                                    
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
    
    const Basis2D&      basis;
    const T             reaction;
    
    Compression1D_t     compression_1d_t;
    Compression1D_x     compression_1d_x;
    Compression2D       compression;
    
    Coefficients<Lexicographical,T,Index2D> P_left_data;
    Coefficients<Lexicographical,T,Index2D> P_right_data;
    
    const LeftPrec2D&   p_left;
    const RightPrec2D&  p_right;
    NoPreconditioner1D  noprec;
    
    const IdentityOperator_t    op_identity_t;
    const IdentityOperator_x    op_identity_x;
    
    const NoInitialCondition    op_noinitcond;
    const InitialCondition&     op_initcond;
    
    DataIdentity_t      data_identity_t;
    DataIdentity_x      data_identity_x;
};    
      
} // namespace lawa

#include <lawa/methods/adaptive/operators/spacetimeoperators/adaptivespacetimereactionoperator1d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVESPACETIMEREACTIONOPERATOR1D_H
 
 

