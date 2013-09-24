/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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

#include <cassert>
#include <lawa/flensforlawa.h>
#include <lawa/math/polynomial.h>
#include <lawa/settings/param.h>
#include <lawa/constructions/realline/subdivision.h>

namespace lawa {

using namespace flens;

template <typename T>
BSpline<T,Dual,Periodic,CDF>::BSpline(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), phiR_(_d, _d_)
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,Periodic,CDF>::BSpline(const MRA<T,Dual,Periodic,CDF> &mra)
    : d(mra.d), d_(mra.d_), mu(d&1), phiR_(d,d_)
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Dual,Periodic,CDF>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(deriv==0);

    if((x < 0.) || (x > 1.)){
        return 0.;
    }

    T val = 0;
    for(int l = ifloor(phiR_.support(j,k).l1); l < iceil<int>(phiR_.support(j,k).l2); ++l){
        val += phiR_(l+x, j, k);
    }
    return val;
}


template <typename T>
PeriodicSupport<T>
BSpline<T,Dual,Periodic,CDF>::support(int j, long k) const
{    
    Support<T> suppR = phiR_.support(j,k);
    if(suppR.length() >= 1){
        return PeriodicSupport<T>(0,1);
    }
    if(suppR.l1 < 0){
        return PeriodicSupport<T>(0,1,suppR.l2, suppR.l1 + 1);
    }
    if(suppR.l2 > 1){
        return PeriodicSupport<T>(0,1,suppR.l2 - 1, suppR.l1);
    }
    return PeriodicSupport<T>(suppR.l1, suppR.l2);
}

template <typename T>
const DenseVector<Array<T> > &
BSpline<T,Dual,Periodic,CDF>::mask() const
{
    return phiR_.a_;
}

} // namespace lawa

