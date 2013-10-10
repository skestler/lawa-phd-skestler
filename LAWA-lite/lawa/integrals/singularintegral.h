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

#ifndef LAWA_INTEGRALS_SINGULARINTEGRAL_H
#define LAWA_INTEGRALS_SINGULARINTEGRAL_H 1

#include <lawa/settings/enum.h>
#include <lawa/functiontypes/function.h>
#include <lawa/integrals/singularquadrature.h>

namespace lawa {

template <typename SingularKernel, typename First, typename Second>
struct SingularIntegral
{
    typedef typename First::T T;

    SingularIntegral(const SingularKernel &singularkernel, const First &first, const Second &second,
                     const T _left=0., const T _right=1.);

    T
    operator()(int _j1, long _k1, XType _e1, int _deriv1,
               int _j2, long _k2, XType _e2, int _deriv2) const;

    T p1(T x) const;

    T p2(T x) const;

    T kernel(T x) const;

    const SingularKernel &singularkernel;
    const First &first;
    const Second &second;
    const T left, right;
    const T RightmLeft;
    const T SqrtRightmLeft;
    const T OneDivSqrtRightmLeft;
    mutable SingularQuadrature<SingularIntegral<SingularKernel,First,Second> > singularquadrature;
    mutable int j1, deriv1,
                j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

template <typename _T>
struct Poly {

    typedef _T T;
    Poly(int _d) : d(_d) {  };

    T
    operator()(T x) const { return std::pow(x,T(d-1)); }

    int d;
};

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
struct SingularIntegralPP
{
    typedef typename FirstPolynomial::T T;

    SingularIntegralPP(const SingularKernel &singularkernel, const FirstPolynomial &first,
                       const SecondPolynomial &second);

    T
    operator()(T a1, T b1, T a2, T b2) const;

    T p1(T x) const;

    T p2(T x) const;

    T kernel(T x) const;

    const SingularKernel   &singularkernel;
    const FirstPolynomial  &first;
    const SecondPolynomial &second;
    mutable SingularQuadrature<SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial> > singularquadrature;
};

} // namespace lawa

#include <lawa/integrals/singularintegral.tcc>

#endif // LAWA_INTEGRALS_SINGULARINTEGRAL_H

