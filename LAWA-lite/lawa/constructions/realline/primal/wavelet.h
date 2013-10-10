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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

using namespace flens;

template <typename _T>
struct Wavelet<_T,Primal,R,CDF>
    : public BasisFunction<_T,Primal,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    Wavelet(int _d, int _d_);

    Wavelet(const BSpline<T,Primal,R,CDF> &_phi,
            const BSpline<T,Dual,R,CDF> &_phi_);

    Wavelet(const Basis<T,Primal,R,CDF> &_basis);

    T
    operator()(T x, int j, long k, unsigned short deriv) const;

    Support<T>
    support(int j, long k) const;

    DenseVector<Array<T> >
    singularSupport(int j, long k) const;

    DenseVector<Array<T> >
    optim_singularSupport(int j, long k) const;

    T
    tic(int j) const;

    const DenseVector<Array<T> > &
    mask() const;

    static DenseVector<Array<T> >
    mask(int d, int d_);

    const int d, d_, mu;
    const int l1, l2;
    const int vanishingMoments;
    const DenseVector<Array<T> > b;
    const BSpline<T,Primal,R,CDF> phi;
    const BSpline<T,Dual,R,CDF> phi_;
    flens::DenseVector<flens::Array<T> > singularPts;
};

} // namespace lawa

#include <lawa/constructions/realline/primal/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_H

