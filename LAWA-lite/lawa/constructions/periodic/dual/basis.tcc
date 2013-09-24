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

namespace lawa {

template <typename T>
Basis<T,Dual,Periodic,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d,d_,j), mra_(d,d_,j), 
      psi_(d,d_), M1_(psi_), _j(j)
{
}

template <typename T>
int
Basis<T,Dual,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Dual,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Dual,Periodic,CDF> &
Basis<T,Dual,Periodic,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi_;
    } else {
        return psi_;
    }
}

template <typename T>
int
Basis<T,Dual,Periodic,CDF>::cardJ_(int j) const
{
    assert(j>=j0);

    return pow2i<T>(j);
}

template <typename T>
Range<int>
Basis<T,Dual,Periodic,CDF>::rangeJ_(int j) const
{
    assert(j>=j0);

    return Range<int>(0,pow2i<T>(j)-1);
}

} // namespace lawa

