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

#include <lawa/aux/compiletime_assert.h>

namespace lawa {

template <FunctionSide Side, Construction Cons, typename X>
typename X::ElementType
evaluate(const MRA<typename X::ElementType,Side,Interval,Cons> &mra, int j,
         const DenseVector<X>& coeffs, typename X::ElementType x, int deriv)
{
    ct_assert(Side==Primal or Side==Orthogonal);

    assert(j>=mra.j0);
    assert(coeffs.length()==mra.cardI(j));
    assert(x>=0.);
    assert(x<=1.);

    typedef typename X::ElementType T;

    BSpline<T,Side,Interval,Cons> phi(mra);
    T ret = 0.0;
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        ret += coeffs(k) * phi(x,j,k,deriv);
    }
    return ret;
/*    if (x<mra.suppIL(j).l2) {
        for (int k=mra.rangeIL(j).firstIndex(); k<=mra.rangeIL(j).lastIndex(); ++k) {
            ret += coeffs(k) * phi(x,j,k);
        }
    }
    if (x>mra.suppIR(j).l1) {
        for (int k=mra.rangeIR(j).firstIndex(); k<=mra.rangeIR(j).lastIndex(); ++k) {
            ret += coeffs(k) * phi(x,j,k);
        }
    }

    int shift = ifloor(pow2i<T>(j) * fabs(mra.phi.support(j,mra.rangeII(j).firstIndex()).l1-x));
    int firstk = mra.rangeII(j).firstIndex();
    int sectork = firstk + shift;

    int start = std::max(sectork-mra.d+1, firstk);
    int end = std::min(sectork+mra.d-1, mra.rangeII(j).lastIndex());
    for (int k=start; k<=end; ++k) {
        ret += coeffs(k) * phi(x,j,k);
    }
    return ret;
*/
}

template <FunctionSide Side, Construction Cons, typename X>
typename X::ElementType
evaluate(const Basis<typename X::ElementType,Side,Interval,Cons> &basis,
         int J, const DenseVector<X> &coeffs, typename X::ElementType x,
         int deriv)
{
    ct_assert(Side==Primal or Side==Orthogonal);

    assert(J>=basis.j0);
    assert(coeffs.range()==basis.mra.rangeI(J));
    assert(x>=0.);
    assert(x<=1.);

    typedef typename X::ElementType T;
/*
    const int j0 = basis.j0;

    basis.setLevel(j0);
    T ret = 0.;
    ret += evaluate(basis.mra,j0,coeffs(basis.mra.rangeI(j0)),x,deriv);
    Wavelet<T,Side,Interval,Cons> psi(basis);
    for (int j=j0; j<=J-1; ++j) {
        if (x<basis.suppJL(j).l2) {
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex() + k) * psi(x, j, k, deriv);
            }
        }
        if (x>basis.suppJR(j).l1) {
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex() + k) * psi(x, j, k);
            }
        }

        T first = basis.suppJI(j).l1;
        T  last = basis.suppJI(j).l2;

        if ((x>first) && (x<last)) {
            int shift = ifloor(pow2i<T>(j) * (x-first));
            int pos = basis.rangeJI(j).firstIndex() + shift;

            int from = std::max(basis.rangeJI(j).firstIndex(),pos-(basis.d+basis.d_));
            int to   = std::min(basis.rangeJI(j).lastIndex(), pos+(basis.d+basis.d_));
            for (int k=from; k<=to; ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex()+k) * psi(x,j,k);
            }
        }
    }
    return ret;
*/
    const int j0 = basis.j0;
    basis.setLevel(j0);

    T ret = 0;
    ret += evaluate(basis.mra,j0,coeffs(Range<int>((int)basis.mra.rangeI(j0).firstIndex(),(int)basis.mra.rangeI(j0).lastIndex())),x,deriv);
    Wavelet<T,Side,Interval,Cons> psi(basis);
    for (int j=j0; j<=J-1; ++j) {
        for (int k=1; k<=basis.cardJ(j); ++k) {
            ret += coeffs(basis.mra.rangeI(j).lastIndex() + k) * psi(x, j, k, deriv);
        }
    }
    return ret;
}

} // namespace lawa

