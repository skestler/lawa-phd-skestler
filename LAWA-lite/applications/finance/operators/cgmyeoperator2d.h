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

#ifndef APPLICATIONS_FINANCE_OPERATORS_CGMYEOPERATOR2D_H
#define APPLICATIONS_FINANCE_OPERATORS_CGMYEOPERATOR2D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/operators/deltas.h>
#include <lawa/settings/settings.h>
//#include <applications/finance/kernels/cgmykernel.h>
//#include <applications/finance/operators/financeoperator1d.h>
#include <applications/finance/operators/cgmyeoperator1d.h>
#include <applications/finance/operators/financeoperator2d.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <typename Basis2D>
struct FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>
{
    typedef typename Basis2D::T T;

    typedef typename Basis2D::FirstBasisType  Basis1;
    typedef typename Basis2D::SecondBasisType Basis2;

    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

    FinanceOperator2D(const Basis2D& _basis,
                      const ProcessParameters2D<T,CGMYeUnivariateJump2D> &_processparameters,
                      T _R1_1=0., T _R2_1=1., T _R1_2=0., T _R2_2=1., int order=10,
                      int _internal_compression_level1=-1, int _internal_compression_level2=-1);

    void
    setCompressionLevel(int _internal_compression_level1, int _internal_compression_level2);

    T
    operator()(const Index2D &row, const Index2D &col);

    void
    eval(Coefficients<Lexicographical,T,Index2D> &v, Coefficients<Lexicographical,T,Index2D> &Av,
         const char* evalType);

    void
    performMV(Coefficients<Lexicographical,T,Index2D> &v, Coefficients<Lexicographical,T,Index2D> &Av);

    const Basis2D                                       &basis;
    const ProcessParameters2D<T,CGMYeUnivariateJump2D>  &processparameters;
    T                                                   R1_1, R2_1, R1_2, R2_2;
    int                                                 internal_compression_level1;
    int                                                 internal_compression_level2;

    FinanceOperator1D<T, CGMYe, Basis1>                 cgmyeop1d_1;
    FinanceOperator1D<T, CGMYe, Basis2>                 cgmyeop1d_2;

    Integral<Gauss, Basis1, Basis1>                     integral1;
    Integral<Gauss, Basis2, Basis2>                     integral2;

    Coefficients<Lexicographical,T,Index2D>             current_v;
    std::map<Index2D,int,lt<Lexicographical,Index2D> >  indicesToInt;
    SparseMatrixT                                       stiffnessMatrix;



};

}

#include <applications/finance/operators/cgmyeoperator2d.tcc>

#endif // APPLICATIONS_FINANCE_OPERATORS_CGMYEOPERATOR2D_H
