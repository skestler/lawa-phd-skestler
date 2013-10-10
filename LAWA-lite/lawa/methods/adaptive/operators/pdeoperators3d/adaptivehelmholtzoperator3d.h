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


#ifndef  LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR3D_H
#define  LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR3D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/compressions/compression_pde1d.h>
#include <lawa/methods/adaptive/compressions/compression_pde3d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
//#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d.h>
#include <lawa/operators/spacetimeoperators/spacetimeoperators.h>
#include <lawa/preconditioners/preconditioners.h>

namespace lawa {


template <typename T, typename Basis3D, typename Preconditioner>
struct AdaptiveHelmholtzOperator3D {
};

template <typename T, typename Basis3D, typename Compression, typename Preconditioner>
struct AdaptiveHelmholtzOperator3D<T, Basis3D, HelmholtzOperator3D<T, Basis3D>, Compression, Preconditioner>
{
    typedef typename Basis2D::FirstBasisType  Basis_x;
    typedef typename Basis2D::SecondBasisType Basis_x;
    typedef typename Basis2D::ThirdBasisType  Basis_z;
    
    
    typedef CompressionPDE1D<T, Basis_x>    Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>    Compression1D_y;
    typedef CompressionPDE1D<T, Basis_z>    Compression1D_z;
    typedef CompressionPDE3D<T, Basis3D>    Compression3D;

    typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;
    
    typedef IdentityOperator1D<T, Basis_x>  IdentityOperator_x;
    typedef IdentityOperator1D<T, Basis_y>  IdentityOperator_y;
    typedef IdentityOperator1D<T, Basis_z>  IdentityOperator_z;
    typedef LaplaceOperator1D<T, Basis_x>   LaplaceOperator_x;
    typedef LaplaceOperator1D<T, Basis_y>   LaplaceOperator_y;
    typedef LaplaceOperator1D<T, Basis_z>   LaplaceOperator_z;

    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataIdentity_x;
    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataIdentity_y;
    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_z,
                               Compression1D_z, NoPreconditioner1D>  DataIdentity_z;
    typedef MapMatrixWithZeros<T, Index1D, LaplaceOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataLaplace_x;
    typedef MapMatrixWithZeros<T, Index1D, LaplaceOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataLaplace_y;
    typedef MapMatrixWithZeros<T, Index1D, LaplaceOperator_z,
                               Compression1D_z, NoPreconditioner1D>  DataLaplace_z;                               
                               

    AdaptiveHelmholtzOperator3D(const Basis3D &_basis, T _c, const Preconditioner &_Prec, 
                                T entrybound=0., int NumOfRows=4096, int NumOfCols=2048);

    T
    operator()(const Index3D &row_index, const Index3D &col_index);

    T
    prec(const Index3D &index);

    void
    clear();
    
    const Basis3D              &basis;
    T                          c;
    const Preconditioner       &Prec;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;
    Compression1D_z            compression_1d_z;
    Compression3D              compression;

    const IdentityOperator_x   op_identity_x;
    const IdentityOperator_y   op_identity_y;
    const IdentityOperator_z   op_identity_z;
    const LaplaceOperator_x    op_laplace_x;
    const LaplaceOperator_y    op_laplace_y;
    const LaplaceOperator_z    op_laplace_z;

    T entrybound;
    int NumOfRows, NumOfCols;

    DataIdentity_x   data_identity_x;
    DataIdentity_y   data_identity_y;
    DataIdentity_z   data_identity_z;
    DataLaplace_x    data_laplace_x;
    DataLaplace_y    data_laplace_y;
    DataLaplace_z    data_laplace_z;

    Coefficients<Lexicographical,T,Index3D> P_data;

};


}    //namespace lawa

#include <lawa/methods/adaptive/operators/adaptivehelmholtzoperator3d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR3D_H

