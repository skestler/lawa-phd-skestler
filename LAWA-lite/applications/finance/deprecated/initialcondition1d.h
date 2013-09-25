/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Sebastian Kestler, Kristina Steih, Mario Rometsch, Alexander Stippler.

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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_INITIALCONDITION1D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_INITIALCONDITION1D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/integrals/integrals.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <applications/finance/options/options.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <OptionType1D OType, typename Basis1D>
struct PayoffInitialCondition1D
{
    typedef typename Basis1D::T T;

    PayoffInitialCondition1D(const Option1D<T,OType> &_option, const Basis1D &_basis,
                             const Function<T> &weight, const T _eta=0,
                             const T _left=0., const T _right=1.);

    T
    operator()(XType _e1, int _j1, long _k1, int _deriv=0) const;

    T
    integrand(T x) const;

    const Option1D<T,OType> &option;
    const Basis1D           &basis;
    const Function<T>       &weight;
    const T                 eta;
    const T                 left, right;
    mutable Quadrature<Gauss,PayoffInitialCondition1D<OType,Basis1D> > quadrature;
    const T                 RightmLeft, SqrtRightmLeft;

    mutable int j1, deriv1,
                j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

}   // namespace lawa

#include <applications/finance/initialconditions/initialcondition1d.tcc>

#endif  // APPLICATIONS_FINANCE_INITIALCONDITIONS_INITIALCONDITION1D_H
