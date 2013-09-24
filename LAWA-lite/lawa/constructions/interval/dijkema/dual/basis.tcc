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

#include <lawa/constructions/interval/initial_stable_completion.h>

namespace lawa {

template <typename T>
Basis<T,Dual,Interval,Dijkema>::Basis(int _d, int _d_, int j)
    : mra(_d, j), mra_(_d, _d_, j),
      d(_d), d_(_d_), mu(d&1),
      min_j0(mra_.min_j0), j0(mra_.j0), _bc(2,0), _j(-1), psi_(*this)
{
    GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
    initial_stable_completion(mra,mra_,Mj1,Mj1_);
    const int cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
    M1_ = RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, Mj1_, min_j0, cons_j);
    _j = std::max(min_j0,j);
    setLevel(_j);
}

template <typename T>
int
Basis<T,Dual,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Dual,Interval,Dijkema>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M1_.setLevel(_j);
    mra.setLevel(_j);
    mra_.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Dual,Interval,Dijkema>::enforceBoundaryCondition()
{
    if ((_bc(0)==0) && (_bc(1)==0)) {
        _bc(0) = _bc(1) = 1;
        mra.enforceBoundaryCondition<BC>();
        mra_.enforceBoundaryCondition<BC>();
        GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
        initial_stable_completion(mra,mra_,Mj1,Mj1_);
        const int cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
        M1_ = RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, 
                                                   Mj1_, min_j0, cons_j);
        setLevel(_j);
    }
}

template <typename T>
const BasisFunction<T,Dual,Interval,Dijkema> &
Basis<T,Dual,Interval,Dijkema>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra_.phi_;
    } else {
        return psi_;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
int
Basis<T,Dual,Interval,Dijkema>::cardJ_(int j) const
{
    assert(j>=j0);
    return pow2i<int>(j);
}

template <typename T>
int
Basis<T,Dual,Interval,Dijkema>::cardJ_L(int j) const
{
    assert(j>=j0);
    return M1_.left.length();
}

template <typename T>
int
Basis<T,Dual,Interval,Dijkema>::cardJ_I(int j) const
{
    assert(j>=j0);
    return M1_.numCols() - M1_.left.length() - M1_.right.length();
}

template <typename T>
int
Basis<T,Dual,Interval,Dijkema>::cardJ_R(int /*j*/) const
{
    return M1_.right.length();
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Dual,Interval,Dijkema>::rangeJ_(int j) const
{
    return Range<int>(1,pow2i<int>(j));
}

template <typename T>
const Range<int>
Basis<T,Dual,Interval,Dijkema>::rangeJ_L(int j) const
{
    assert(j>=j0);
    return Range<int>(1,M1_.left.length());
}

template <typename T>
const Range<int>
Basis<T,Dual,Interval,Dijkema>::rangeJ_I(int j) const
{
    assert(j>=j0);
    return Range<int>(M1_.left.length()+1, pow2i<int>(j)-M1_.right.length());
}

template <typename T>
const Range<int>
Basis<T,Dual,Interval,Dijkema>::rangeJ_R(int j) const
{
    assert(j>=j0);
    return Range<int>(pow2i<int>(j)-M1_.right.length()+1, pow2i<int>(j));
}

} // namespace lawa

