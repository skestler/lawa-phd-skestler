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
BSpline<T,Dual,Interval,Dijkema>::BSpline(const MRA<T,Dual,Interval,Dijkema> &_mra_)
    : mra_(_mra_)
{
}

template <typename T>
T
BSpline<T,Dual,Interval,Dijkema>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(j>=mra_.j0);
    assert(k>=mra_.rangeI_(j).firstIndex());
    assert(k<=mra_.rangeI_(j).lastIndex());
    assert(deriv==0);

    if (mra_._bc(0)==1 && mra_._bc(1)==1 && mra_.d!=3 && mra_.d_!=5) {
        std::cerr << "BSpline<T,Dual,Interval,Dijkema>::operator(...): here is still an error in the "
                     "evaluation of the right boundary functions." << std::endl;
    }

    //Todo: There is something wrong with the indexing when Dirichlet bc are enforced (except for
    //      (d,d_)=(3,5) at the right boundary.
    if (k<=mra_.rangeI_L().lastIndex()) {
        T value = 0.0;
        int l = -mra_.phi_R.a_.lastIndex()+1;
        for (int r=mra_.R_Left.firstRow(); r<=mra_.R_Left.lastRow(); ++r, ++l) {
            //std::cerr << "Left: k = " << k << ", l = " << l  << ": " << mra_.R_Left(r,k)<< std::endl;
            value += mra_.R_Left(r,k) * mra_.phi_R(x,j,l);
        }
//        assert(l==-mra_.phi_R.a_.firstIndex());
        return value;
    }

    if (k>=mra_.rangeI_R(j).firstIndex()) {
        k -= mra_.rangeI_R(j).firstIndex()-1;
        T value = 0.0;
        int l = pow2i<T>(j)-mra_.phi_R.a_.lastIndex()+1;
        for (int r=mra_.R_Right.firstRow(); r<=mra_.R_Right.lastRow(); ++r, ++l) {
            //std::cerr << "Right: k = " << k << ", l = " << l  << ": " << mra_.R_Right(r,k) << std::endl;
            value += mra_.R_Right(r,k) * mra_.phi_R(x,j,l);
        }
        return value;
    }

    return mra_.phi_R(x,j,k-(mra_.d + mra_.d_ - 1)-mra_.phi_R.l1_);
}

template <typename T>
Support<T>
BSpline<T,Dual,Interval,Dijkema>::support(int j, long k) const
{
    if (k<=mra_.rangeI_L().lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j)*(mra_.phi_R.a_.length()-2));
    }
    if (k>=mra_.rangeI_R(j).firstIndex()) {
        return Support<T>(1.-pow2i<T>(-j)*(mra_.phi_R.a_.length()-2), 1.);
    }
    return pow2i<T>(-j)*Support<T>(k-(mra_.d+mra_.d_-1),
                                   k-(mra_.d+mra_.d_-1)+mra_.phi_R.l2_-mra_.phi_R.l1_);
}

} // namespace lawa

