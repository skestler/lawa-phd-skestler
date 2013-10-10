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


#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEIDENTITYOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEIDENTITYOPERATOR2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>

namespace lawa {

template <typename T, typename Basis2D>
struct AdaptiveIdentityOperator2D
{
    typedef typename Basis2D::FirstBasisType  Basis_x;
    typedef typename Basis2D::SecondBasisType Basis_y;

    typedef NoPreconditioner<T,Index1D>                              NoPreconditioner1D;
    typedef NoPreconditioner<T,Index2D>                              NoPreconditioner2D;

    typedef CompressionPDE1D<T, Basis_x>                             Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>                             Compression1D_y;
    typedef CompressionPDE2D<T, Basis2D>                             Compression2D;

    typedef IdentityOperator1D<T, Basis_x>                           IdentityOperator_x;
    typedef IdentityOperator1D<T, Basis_y>                           IdentityOperator_y;

    typedef MapMatrix<T, Index1D, IdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataReaction_x;
    typedef MapMatrix<T, Index1D, IdentityOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataReaction_y;

    AdaptiveIdentityOperator2D(const Basis2D &_basis, T _entrybound=0.,
                                   int _NumOfRows=4096, int _NumOfCols=2048);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    const Basis2D                     &basis;

    Compression1D_x                   compression_1d_x;
    Compression1D_y                   compression_1d_y;
    Compression2D                     compression;
    const NoPreconditioner2D          Prec;
    const IdentityOperator_x          op_identity_x;
    const IdentityOperator_y          op_identity_y;

    T entrybound;
    int NumOfRows, NumOfCols;

    DataReaction_x   data_identity_x;
    DataReaction_y   data_identity_y;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptiveidentityoperator2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEIDENTITYOPERATOR2D_H

