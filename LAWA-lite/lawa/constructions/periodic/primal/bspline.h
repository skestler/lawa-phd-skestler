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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/periodic/periodicsupport.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Primal,Periodic,CDF>
    : public BasisFunction<_T,Primal,Periodic,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = Periodic;
    static const Construction Cons = CDF;

    BSpline(int _d);                        // Periodized Bspline on [0,1]

    BSpline(const MRA<T,Primal,Periodic,CDF> &mra);

    virtual
    ~BSpline();

    T
    operator()(T x, int j, long k, unsigned short deriv) const;    // (cumulated) evaluation on [0,1]

    PeriodicSupport<T>
    support(int j, long k) const;            // subset of [0,1], may be divided in 2 parts

    DenseVector<Array<T> >
    singularSupport(int j, long k) const;    // Points in [0,1]

    T
/**/tic(int j) const;                       // -> ???

    const DenseVector<Array<T> > &
    mask() const;                           // mask of original Spline on R

    const int d, mu;
    const BSpline<T, Primal, R, CDF> phiR;  // ''original'' Bspline on R    
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H

