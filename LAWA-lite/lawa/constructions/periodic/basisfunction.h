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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_BASISFUNCTION_H
#define LAWA_CONSTRUCTIONS_PERIODIC_BASISFUNCTION_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/periodic/periodicsupport.h>

namespace lawa {

template <typename _T, FunctionSide _Side, Construction _Cons>
struct BasisFunction<_T, _Side, Periodic, _Cons>
{
    typedef _T T;
    static const FunctionSide Side = _Side;
    static const DomainType Domain = Periodic;
    static const Construction Cons = _Cons;

    virtual T
    operator()(T x, int j, long k, unsigned short deriv) const;

    virtual PeriodicSupport<T>
    support(int j, long k) const;

    virtual DenseVector<Array<T> >
    singularSupport(int j, long k) const;

    virtual T
    tic(int j) const;
};

} // namespace lawa

#include <lawa/constructions/periodic/basisfunction.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_BASISFUNCTION_H

