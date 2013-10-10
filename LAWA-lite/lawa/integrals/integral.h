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

#ifndef LAWA_INTEGRALS_INTEGRAL_H
#define LAWA_INTEGRALS_INTEGRAL_H 1

#include <lawa/settings/enum.h>
#include <lawa/functiontypes/function.h>
#include <lawa/integrals/quadrature.h>

namespace lawa {

template <QuadratureType Quad, typename First, typename Second>
struct Integral
{
    typedef typename First::T T;

    Integral(const First &first, const Second &second);

    T
    operator()(int _j1, long _k1, XType _e1, int _deriv1,
               int _j2, long _k2, XType _e2, int _deriv2) const;

    T
    integrand(T x) const;

    const First &first;
    const Second &second;
    T left, right;
    mutable Quadrature<Quad,Integral<Quad,First,Second> > quadrature;
    mutable int j1, deriv1,
                j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

template <QuadratureType Quad, typename First, typename Second=First>
struct IntegralF
{
    typedef typename First::T T;

    IntegralF(const Function<T> &_function, const First &_first,
              const T _left=0., const T _right=1.);

    IntegralF(const Function<T> &_function, const First &_first, const Second &_second,
              const T _left=0., const T _right=1.);

    T
    operator()(int _j1, long _k1, XType _e1, int _deriv) const;

    T
    operator()(int _j1, long _k1, XType _e1, int _deriv1,
               int _j2, long _k2, XType _e2, int _deriv2) const;

    T
    integrand(T x) const;

    const Function<T> &function;
    const First &first;
    const Second &second;
    const T left;
    const T right;
    const T RightmLeft;
    const T SqrtRightmLeft;
    const T OneDivSqrtRightmLeft;
    Quadrature<Quad,IntegralF<Quad,First,Second> > quadrature;
    const bool _f2;
    mutable int j1, deriv1,
                j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

} // namespace lawa

#include <lawa/integrals/integral.tcc>

#endif // LAWA_INTEGRALS_INTEGRAL_H

