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

#ifndef LAWA_INTEGRALS_INTEGRAL2D_H
#define LAWA_INTEGRALS_INTEGRAL2D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/functiontypes/function2d.h>
#include <lawa/integrals/quadrature2d.h>

namespace lawa {

template <QuadratureType Quad, typename BasisX, typename BasisY>
struct Integral2D
{
    typedef typename BasisX::T T;
    typedef BasisX Basis_X;
    typedef BasisY Basis_Y;
    
    const Function2D<T> &F;
    const BasisX &basisx;
    const BasisY &basisy;
    mutable int jx, jy;
    mutable unsigned short derivx, derivy;
    mutable long kx, ky;
    mutable XType ex, ey;
    Quadrature2D<Quad, Integral2D<Quad,BasisX,BasisY> > quadrature;

    Integral2D(const Function2D<T> &_F, const BasisX &basisx, const BasisY &basisy);

    T
    operator()(int _jx, long _kx, XType _ex, unsigned short _derivx,
               int _jy, long _ky, XType _ey, unsigned short _derivy) const;

    T
    integrand(T x, T y) const;
};

} // namespace lawa

#include <lawa/integrals/integral2d.tcc>

#endif  //LAWA_INTEGRALS_INTEGRAL2D_H

