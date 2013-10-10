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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_EUROPEANBASKET2D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_EUROPEANBASKET2D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/integrals/integrals.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <applications/finance/initialconditions/adaptivepayoffquadrature2d.h>
#include <applications/finance/options/options.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
struct PayoffIntegral2D
{
    typedef typename Basis::T T;
    typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

    PayoffIntegral2D(const Basis &_basis, const PayoffFunction &_payofffunction,
                     const T _left_x1=0., const T _right_x1=1.,
                     const T _left_x2=0., const T _right_x2=1.,
                     bool _useSpecialRefinement=false, T _maxRectangleLength=0.25, int _order=20);

    T
    integrand(T x1, T x2) const;

    T
    operator()(const Index2D &lambda) const;

    T
    integrate(T a1, T b1, T a2, T b2) const;

    T
    _getOrderAndValue(T a1, T b1, T a2, T b2) const;

    const Basis                     &basis;
    const PayoffFunction            &payofffunction;

    Quadrature2D<Quad, PayoffIntegral2D<Quad,Basis,PayoffFunction> > highestOrderQuadrature;
    Quadrature2D<Quad, PayoffIntegral2D<Quad,Basis,PayoffFunction> > highOrderQuadrature;
    Quadrature2D<Quad, PayoffIntegral2D<Quad,Basis,PayoffFunction> > mediumOrderQuadrature;
    Quadrature2D<Quad, PayoffIntegral2D<Quad,Basis,PayoffFunction> > lowOrderQuadrature;

    const T                         left_x1, right_x1;
    const T                         left_x2, right_x2;
    const T                         RightmLeft_x1, SqrtRightmLeft_x1;
    const T                         RightmLeft_x2, SqrtRightmLeft_x2;

    bool                            useSpecialRefinement;
    T                               maxRectangleLength;
    int                             order;

    mutable int j1, deriv1, j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

}   // namespace lawa

#include <applications/finance/initialconditions/payoffintegral2d.tcc>

#endif  // APPLICATIONS_FINANCE_INITIALCONDITIONS_BASKETPUT2D_H
