/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

namespace lawa {
    

template <typename X>
typename X::ElementType
evaluate(const MRA<typename X::ElementType,Primal,Periodic,CDF> &mra, int j,
         const DenseVector<X> &coeffs, typename X::ElementType x, int deriv)
{
    typedef typename X::ElementType T;
    assert(j>=mra.j0);
    assert(coeffs.length()==mra.cardI(j));
    assert(x>=0.);
    assert(x<=1.);
    
    BSpline<T,Primal,Periodic,CDF> phi(mra.d);
    int offsetI = mra.rangeI(mra.j0).firstIndex()-coeffs.firstIndex();
    T ret = 0.0;
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        ret += coeffs(k-offsetI) * phi(x,j,k,deriv);
    }
    return ret;
}

template <typename X>
typename X::ElementType
evaluate(const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis,
         int J, const DenseVector<X> &coeffs, typename X::ElementType x, 
         int deriv)
{
    typedef typename X::ElementType T;
    assert(J>=basis.j0);
    assert(coeffs.length()==basis.mra.cardI(J));
    assert(x>=0.);
    assert(x<=1.);

    const int j0 = basis.j0;
    basis.setLevel(j0);
    T ret = 0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-coeffs.firstIndex();
    Range<int> range(coeffs.firstIndex(), coeffs.firstIndex() + basis.cardJ(j0) - 1);
    
    ret += evaluate(basis.mra,j0,coeffs(range),x,deriv);
    Wavelet<T,Primal,Periodic,CDF> psi(basis.d, basis.d_);
    for (int j=j0; j<=J-1; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            ret += coeffs(basis.mra.cardI(j) + k - offsetJ) * basis.psi(x,j,k,deriv);
        }
    }
    return ret;
}

} // namespace lawa

