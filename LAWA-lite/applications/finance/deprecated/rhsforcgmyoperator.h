/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef APPLICATIONS_FINANCE_RHSFORCGMYOPERATOR_H
#define APPLICATIONS_FINANCE_RHSFORCGMYOPERATOR_H 1

#include <applications/finance/cgmyutils.h>

namespace lawa {

template <typename T, typename Basis>
class RHSForCGMYOperator1D
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

public:
    RHSForCGMYOperator1D(const Basis &basis, T (*f)(T), T (*df)(T),
                         const DenseVector<Array<T> > &_singularPoints,
                         T C, T G, T M, T Y, int _order);

    T
    operator()(XType xtype, int j, int k) const;

    T
    operator()(const Index1D &lambda) const;

private:
    flens::GeMatrix<FullStorage<T,ColMajor> >
    _computeDeltas(XType xtype, int j, int k) const;

    T
    _evaluate_acal_u(T x_ast, T a, T b) const;

    T
    _quadrature_du_vs_CGMYkernel(T x_ast, T a, T b) const;

    T
    _quadrature_ddu_vs_CGMYkernel(T x_ast, T a, T b) const;

    const Basis &basis;
    PrimalSpline phi;
    PrimalWavelet psi;
    Function<T> f, df;
    CGMYUtils<T> cgmy;
    int order;
    flens::GeMatrix<FullStorage<T,ColMajor> > phi_deltas;
    flens::GeMatrix<FullStorage<T,ColMajor> > psi_deltas;
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;

};

}    //namespace lawa

#include <applications/finance/rhsforcgmyoperator.tcc>

#endif    //APPLICATIONS_FINANCE_RHSFORCGMYOPERATOR_H
