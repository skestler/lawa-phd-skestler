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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/periodic/periodicsupport.h>
#include <lawa/constructions/periodic/primal/bspline.h>
#include <lawa/constructions/periodic/dual/bspline.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

using namespace flens;

template <typename _T>
struct Wavelet<_T,Dual,Periodic,CDF>
    : public BasisFunction<_T,Dual,Periodic,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Dual;
    static const DomainType Domain = Periodic;
    static const Construction Cons = CDF;

    Wavelet(int _d, int _d_);

    Wavelet(const BSpline<T,Primal,Periodic,CDF> &_phi,
            const BSpline<T,Dual,Periodic,CDF> &_phi_);

    Wavelet(const Basis<T,Dual,Periodic,CDF> &_basis);

    T
    operator()(T x, int j, long k, unsigned short deriv=0) const;

    PeriodicSupport<T>
    support(int j, long k) const;

    const DenseVector<Array<T> > &
    mask() const;

    static DenseVector<Array<T> >
    mask(int d, int d_);

    const int d, d_, mu;
    const Wavelet<T, Dual, R, CDF> psiR_;
};

} // namespace lawa

#include <lawa/constructions/periodic/dual/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_H

