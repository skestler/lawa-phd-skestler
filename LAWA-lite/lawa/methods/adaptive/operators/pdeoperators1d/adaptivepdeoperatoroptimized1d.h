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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEPDEOPERATOROPTIMIZED1D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEPDEOPERATOROPTIMIZED1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/pdeoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct AdaptivePDEOperatorOptimized1D {

};


template <typename T>
struct AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF> {

    typedef Basis<T,Primal,R,CDF>                                             ReallineCDFBasis1D;

    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;

    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::iterator        coeff1d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator    const_abs_coeff1d_it;

    typedef CompressionPDE1D<T, ReallineCDFBasis1D>                           Compression1D;

    typedef NoPreconditioner<T,Index1D>                                       NoPreconditioner1D;

    typedef PDEOperator1D<T, ReallineCDFBasis1D>                              PDEOp1D;

    typedef MapMatrix<T, Index1D, PDEOp1D,
                     Compression1D, NoPreconditioner1D>                       DataPDE1D;

    AdaptivePDEOperatorOptimized1D(const ReallineCDFBasis1D &_basis, T _reaction, T _convection,
                                   T _diffusion, T thresh=0., int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    prec(const Index1D &index);

    Coefficients<Lexicographical,T,Index1D>
    mv(const IndexSet<Index1D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index1D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, T eps);


    Coefficients<Lexicographical,T,Index1D>
    apply(const Coefficients<Lexicographical,T,Index1D> &v, int k=0, int J=-1000);

    void
    apply(const Coefficients<Lexicographical,T,Index1D> &v, T eps,
          Coefficients<Lexicographical,T,Index1D> &ret);

    void
    apply(const Coefficients<Lexicographical,T,Index1D> &v, T eps,
          const IndexSet<Index1D> &Lambda, Coefficients<Lexicographical,T,Index1D> &ret);

    int
    findK(const Coefficients<AbsoluteValue,T,Index1D> &v, T eps);

    void
    clear();

    const ReallineCDFBasis1D                &basis;
    T                                       reaction, convection, diffusion;
    T cA, CA, kappa;

    Compression1D                           compression;
    const PDEOp1D                           pde_op1d;
    NoPreconditioner1D                      prec1d;
    DataPDE1D                               pde_data1d;

    Coefficients<Lexicographical,T,Index1D> P_data;
};



template <typename T, DomainType Domain>
struct AdaptivePDEOperatorOptimized1D<T,Primal,Domain,SparseMulti>
{
    typedef Basis<T,Primal,Domain,SparseMulti>  SparseMultiBasis1D;
    ct_assert( IsRealline<SparseMultiBasis1D>::value);


    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::iterator        coeff1d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator    const_abs_coeff1d_it;

    typedef CompressionPDE1D<T, SparseMultiBasis1D>                           Compression1D;

    typedef NoPreconditioner<T,Index1D>                                       NoPreconditioner1D;

    typedef PDEOperator1D<T, SparseMultiBasis1D>                              PDEOp1D;

    typedef MapMatrix<T, Index1D, PDEOp1D,
                      Compression1D, NoPreconditioner1D>                      DataPDE1D;

    AdaptivePDEOperatorOptimized1D(const SparseMultiBasis1D &_basis1d, T _reaction, T _convection,
                                   T _diffusion, T thresh=0., int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    prec(const Index1D &index);

    Coefficients<Lexicographical,T,Index1D>
    mv(const IndexSet<Index1D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index1D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, T eps);

    Coefficients<Lexicographical,T,Index1D>
    apply(const Coefficients<Lexicographical,T,Index1D> &v, int k=0, int J=0,
          cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    apply(const Coefficients<Lexicographical,T,Index1D> &v, T eps,
          Coefficients<Lexicographical,T,Index1D> &ret, cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    apply(const Coefficients<Lexicographical,T,Index1D> &v, T eps,
          const IndexSet<Index1D> &Lambda, Coefficients<Lexicographical,T,Index1D> &ret,
          cxxblas::Transpose trans=cxxblas::NoTrans);

    void
    clear();

    const SparseMultiBasis1D     &basis;
    T reaction, convection, diffusion;
    T cA, CA, kappa;

    Compression1D                compression;
    const PDEOp1D                pde_op1d;
    NoPreconditioner1D           prec1d;
    DataPDE1D                    pde_data1d;

    Coefficients<Lexicographical,T,Index1D> P_data;
};

}

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivepdeoperatoroptimized1d.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEPDEOPERATOROPTIMIZED1D_H
