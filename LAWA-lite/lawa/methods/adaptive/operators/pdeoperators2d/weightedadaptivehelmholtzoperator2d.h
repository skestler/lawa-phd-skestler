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


#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_WEIGHTEDADAPTIVEHELMHOLTZOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_WEIGHTEDADAPTIVEHELMHOLTZOPERATOR2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/compressions/compression_weightedpde1d.h>
#include <lawa/methods/adaptive/compressions/compression_weightedpde2d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/weightedidentityoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedlaplaceoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>
#include <lawa/preconditioners/preconditioners.h>

namespace lawa {

template <typename T, typename Basis, typename Preconditioner>
struct WeightedAdaptiveHelmholtzOperator2D : public Operator2D<T>
{
    typedef typename Basis::FirstBasisType  Basis_x;
    typedef typename Basis::SecondBasisType Basis_y;

    typedef CompressionWeightedPDE1D<T, Basis_x>            Compression1D_x;
    typedef CompressionWeightedPDE1D<T, Basis_y>            Compression1D_y;
    typedef CompressionWeightedPDE2D<T, Basis>         	    Compression2D;

    typedef NoPreconditioner<T,Index1D>             NoPreconditioner1D;
    typedef NoPreconditioner<T,Index2D>             NoPreconditioner2D;

    typedef WeightedIdentityOperator1D<T, Basis_x, Gauss>   WeightedIdentityOperator_x;
    typedef WeightedIdentityOperator1D<T, Basis_y, Gauss>   WeightedIdentityOperator_y;
    typedef WeightedLaplaceOperator1D<T, Basis_x, Gauss>    WeightedLaplaceOperator_x;
    typedef WeightedLaplaceOperator1D<T, Basis_y, Gauss>    WeightedLaplaceOperator_y;

    typedef MapMatrix<T, Index1D, WeightedIdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataWeightedIdentity_x;
    typedef MapMatrix<T, Index1D, WeightedIdentityOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataWeightedIdentity_y;
    typedef MapMatrix<T, Index1D, WeightedLaplaceOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataWeightedLaplace_x;
    typedef MapMatrix<T, Index1D, WeightedLaplaceOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataWeightedLaplace_y;

    WeightedAdaptiveHelmholtzOperator2D(const Basis &_basis, const T _c, 
                                      Function<T> weightFct_x, Function<T> weightFct_y,
                                      const Preconditioner &_Prec);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    void
    clear();


    const Basis&            basis;
    const T                 c;
    const Preconditioner&   Prec;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;
    Compression2D              compression;

    const WeightedIdentityOperator_x    op_identity_x;
    const WeightedIdentityOperator_y    op_identity_y;
    const WeightedLaplaceOperator_x     op_laplace_x;
    const WeightedLaplaceOperator_y     op_laplace_y;

    NoPreconditioner1D         Prec1D;
    
    DataWeightedIdentity_x    data_identity_x;
    DataWeightedIdentity_y    data_identity_y;
    DataWeightedLaplace_x     data_laplace_x;
    DataWeightedLaplace_y     data_laplace_y;

    Coefficients<Lexicographical,T,Index2D> P_data;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators2d/weightedadaptivehelmholtzoperator2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_WEIGHTEDADAPTIVELHELMHOLTZOPERATOR2D_H

