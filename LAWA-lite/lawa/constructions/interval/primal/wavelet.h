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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_WAVELET_H 1

#include <lawa/constructions/support.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T, Construction _Cons>
struct Wavelet<_T,Primal,Interval,_Cons>
    : public BasisFunction<_T,Primal,Interval,_Cons>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = Interval;
    static const Construction Cons = _Cons;

    Wavelet(const Basis<T,Primal,Interval,Cons> &_basis);

    T
    operator()(T x, int j, long k, unsigned short deriv) const;

    Support<T>
    support(int j, long k) const;

    DenseVector<Array<T> >
    singularSupport(int j, long k) const;

    int
    vanishingMoments(int j, long k) const;

    T
    tic(int j) const;

    DenseVector<Array<long double> > *
    getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const;

    int
    getRefinementLevel(int j) const;

    // returns the range of B-splines required for the refinement of $\psi_{j,k}$
    void
    getRefinementNeighbors(int j, long k, int &refinement_j, long &refinement_k_first,
                           long &refinement_k_last) const;

    T
    getL2Norm(int j, long k) const;

    T
    getH1SemiNorm(int j, long k) const;

    const Basis<T,Primal,Interval,Cons> &basis;
};

} // namespace lawa

#include <lawa/constructions/interval/primal/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_WAVELET_H

