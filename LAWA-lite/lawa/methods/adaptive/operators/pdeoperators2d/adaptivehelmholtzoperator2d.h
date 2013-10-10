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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptiveoperator2d.h>

namespace lawa {

template <typename T, typename Basis2D, typename Preconditioner>
struct AdaptiveHelmholtzOperator2D : public AdaptiveOperator2D<T>
{
    typedef typename Basis2D::FirstBasisType  Basis_x;
    typedef typename Basis2D::SecondBasisType Basis_y;

    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > SparseMatrixT;

    typedef CompressionPDE1D<T, Basis_x>                             Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>                             Compression1D_y;
    typedef CompressionPDE2D<T, Basis2D>                             Compression2D;

    typedef NoPreconditioner<T,Index1D>                              NoPreconditioner1D;
    typedef NoPreconditioner<T,Index2D>                              NoPreconditioner2D;

    typedef IdentityOperator1D<T, Basis_x>                           IdentityOperator_x;
    typedef IdentityOperator1D<T, Basis_y>                           IdentityOperator_y;
    typedef LaplaceOperator1D<T, Basis_x>                            LaplaceOperator_x;
    typedef LaplaceOperator1D<T, Basis_y>                            LaplaceOperator_y;

    typedef MapMatrix<T, Index1D, IdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataIdentity_x;
    typedef MapMatrix<T, Index1D, IdentityOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataIdentity_y;
    typedef MapMatrix<T, Index1D, LaplaceOperator_x,
                               Compression1D_x, NoPreconditioner1D>  DataLaplace_x;
    typedef MapMatrix<T, Index1D, LaplaceOperator_y,
                               Compression1D_y, NoPreconditioner1D>  DataLaplace_y;

    typedef HelmholtzOperator2D<T, Basis2D>                 HelmholtzBilinearForm2D;
    typedef DiagonalMatrixPreconditioner2D<T, Basis2D,
                                        HelmholtzBilinearForm2D>  DiagonalHelmholtzPreconditioner2D;

    AdaptiveHelmholtzOperator2D(const Basis2D &_basis2d, T _c, const Preconditioner &_Prec,
                                T thresh=0., int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1, bool useLinearIndex=false);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, T eps, bool useLinearIndex=false);

    Coefficients<Lexicographical,T,Index2D>
    apply(const Coefficients<Lexicographical,T,Index2D> &v, int k, int J=-1000,
          cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
          Coefficients<Lexicographical,T,Index2D> &ret,
          cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
          const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
          cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    clear();


    const Basis2D              &basis;
    T c;
    const Preconditioner       &Prec;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;
    Compression2D              compression;

    const IdentityOperator_x   op_identity_x;
    const IdentityOperator_y   op_identity_y;
    const LaplaceOperator_x    op_laplace_x;
    const LaplaceOperator_y    op_laplace_y;

    NoPreconditioner1D         Prec1D;


    DataIdentity_x   data_identity_x;
    DataIdentity_y   data_identity_y;
    DataLaplace_x    data_laplace_x;
    DataLaplace_y    data_laplace_y;

    Coefficients<Lexicographical,T,Index2D> P_data;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptivehelmholtzoperator2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOR2D_H

