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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOROPTIMIZED2D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOROPTIMIZED2D_H 1

#include <ext/hash_map>

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/preconditioners/preconditioners.h>
//#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivelaplaceoperator1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveidentityoperator1d.h>

namespace lawa {

template <typename T, FunctionSide Side1, DomainType Domain1, Construction Cons1,
                      FunctionSide Side2, DomainType Domain2, Construction Cons2>
struct AdaptiveHelmholtzOperatorOptimized2D
{

};

template <typename T, DomainType Domain>
struct AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
{
    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef IndexSet<Index2D>::const_iterator                                 const_set2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator  const_coeff2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::iterator        coeff2d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index2D>::const_iterator    const_abs_coeff2d_it;

    typedef Basis<T,Orthogonal,Domain,Multi>                                 Basis_x;
    typedef Basis<T,Orthogonal,Domain,Multi>                                 Basis_y;
    typedef TensorBasis2D<Adaptive,Basis_x,Basis_y>                           Basis2D;

    //ct_assert(   IsRealline<typename Basis2D::FirstBasisType>::value
    //          && IsRealline<typename Basis2D::SecondBasisType>::value);

    typedef CompressionPDE1D<T, Basis_x>                                      Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>                                      Compression1D_y;
    typedef CompressionPDE2D<T, Basis2D>                                      Compression2D;

    typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Domain,Multi>              DataLaplace1D;

    AdaptiveHelmholtzOperatorOptimized2D(const Basis2D &_basis2d, T _c);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

    /*
    Coefficients<Lexicographical,T,Index2D>
    operator*(const Coefficients<Lexicographical,T,Index2D> &v);
    */
    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, T eps);

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

    int
    findK(const Coefficients<AbsoluteValue,T,Index2D> &v, T eps);

    void
    clear();

    const Basis2D              &basis;
    T c;

    T cA, CA, kappa;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;
    Compression2D              compression;

    DataLaplace1D    laplace_data1d;

    Coefficients<Lexicographical,T,Index2D> P_data;
};

template <typename T, DomainType Domain1, DomainType Domain2>
struct AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
{
    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef IndexSet<Index2D>::const_iterator                                 const_set2d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::iterator        coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator  const_coeff2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::iterator        coeff2d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index2D>::const_iterator    const_abs_coeff2d_it;

    /*
    typedef typename __gnu_cxx::hash_map<Index1D,
                                         Coefficients<Lexicographical,T,Index1D>,
                                         index_hashfunction<Index1D>,
                                         index_eqfunction<Index1D> >          Index1D_Coefficients1D_Hash;
     */
    typedef typename std::map<Index1D, Coefficients<Lexicographical,T,Index1D>,
                              lt<Lexicographical, Index1D> >                       Index1D_Coefficients1D_Hash;
    typedef typename Index1D_Coefficients1D_Hash::iterator                    Index1D_Coefficients1D_Hash_it;
    typedef typename Index1D_Coefficients1D_Hash::const_iterator              const_Index1D_Coefficients1D_Hash_it;

    typedef Basis<T,Primal,Domain1,SparseMulti>                               Basis_x;
    typedef Basis<T,Primal,Domain2,SparseMulti>                               Basis_y;
    typedef TensorBasis2D<Adaptive,Basis_x,Basis_y>                           Basis2D;

    typedef CompressionPDE1D<T, Basis_x>                                      Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>                                      Compression1D_y;
    typedef CompressionPDE2D<T, Basis2D>                                      Compression2D;

    typedef AdaptiveLaplaceOperator1D<T,Primal,Domain1,SparseMulti>           DataLaplace1D_x;
    typedef AdaptiveIdentityOperator1D<T,Primal,Domain1,SparseMulti>          DataIdentity1D_x;
    typedef AdaptiveLaplaceOperator1D<T,Primal,Domain2,SparseMulti>           DataLaplace1D_y;
    typedef AdaptiveIdentityOperator1D<T,Primal,Domain2,SparseMulti>          DataIdentity1D_y;

    AdaptiveHelmholtzOperatorOptimized2D(const Basis2D &_basis2d, T _c, T _diffusion_y=1.);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, T eps);

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
    T diffusion_y;
    T thresh;

    T cA, CA, kappa;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;
    Compression2D              compression;

    DataLaplace1D_x            laplace_data1d_x;
    DataIdentity1D_x           identity_data1d_x;
    DataLaplace1D_y            laplace_data1d_y;
    DataIdentity1D_y           identity_data1d_y;

    Coefficients<Lexicographical,T,Index2D> P_data;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptivehelmholtzoperatoroptimized2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATOROPTIMIZED2D_H

