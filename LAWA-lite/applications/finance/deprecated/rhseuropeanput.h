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

#ifndef APPLICATIONS_FINANCE_RHSEUROPEANPUT_H
#define APPLICATIONS_FINANCE_RHSEUROPEANPUT_H 1

#include <lawa/integrals/integrals.h>

namespace lawa {

template <typename T>
struct StabilizingExpFunction
{
    static T lambda;
    static DenseVector<Array<T> > sing_pts;
    static void setlambda(T _lambda);
    static T ExpPowMlambda(T t);
};

template <typename T, typename Basis>
class RHSEuropeanPut
{
    typedef typename Basis::FirstBasisType Basis_t;
    typedef typename Basis::SecondBasisType Basis_x;

    typedef typename Basis_t::BSplineType PrimalSpline_t;
    typedef typename Basis_x::BSplineType PrimalSpline_x;
    typedef typename Basis_t::WaveletType PrimalWavelet_t;
    typedef typename Basis_x::WaveletType PrimalWavelet_x;

public:
    RHSEuropeanPut(const Basis &basis, T (*f)(T), DenseVector<Array<T> > &sing_pts,
                   T _C, T _G, T _M, T _Y, T _r, T _K, T _maturity, T _sigma=0.);

    T
    operator()(XType xtype, int j, int k) const;

    T
    operator()(XType xtype_t, int j_t, int k_t, XType xtype_x, int j_x, int k_x) const;

    T
    operator()(const Index1D &lambda) const;

    T
    operator()(const Index2D &lambda) const;

private:
    const Basis &basis;
    Function<T> stabilizing_exp_function;

    PrimalSpline_t  phi_t;
    PrimalWavelet_t psi_t;
    PrimalSpline_x  phi_x, d_phi_x;
    PrimalWavelet_x psi_x, d_psi_x;
    Integral<T, Gauss, PrimalSpline_t,  Function<T> > integral_sff;
    Integral<T, Gauss, PrimalWavelet_t, Function<T> > integral_wf;

    CGMYUtils<T> cgmy;
    T r;
    T K;
    T maturity;
    T sigma;

    CGMYOperator1D_Fast<T, Basis_x> cgmy_operator;

};


}    //namespace lawa

#include <applications/finance/rhseuropeanput.tcc>

#endif //APPLICATIONS_FINANCE_RHSEUROPEANPUT_H
