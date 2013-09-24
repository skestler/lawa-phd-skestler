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

#include <cassert>

namespace lawa {

template <typename X, typename Y, Construction Cons>
void
decompose(const DenseVector<X> &x, 
          const Basis<typename X::ElementType,Dual,Interval,Cons> &basis_, int j,
          DenseVector<Y> &y)
{
    assert(j>=basis_.j0);
    assert(x.range()==basis_.mra_.rangeI_(j+1));
    typedef typename X::ElementType T;

    basis_.setLevel(j);
    DenseVector<Array<T> > sf, w;
    sf = transpose(basis_.mra_.M0_) * x;
    w  = transpose(basis_.M1_) * x;
    concatenate(sf,w,y,y.firstIndex());
}

template <typename X, typename Y, Construction Cons>
void
decompose_(const DenseVector<X> &x,
           const Basis<typename X::ElementType,Primal,Interval,Cons> &basis, int j,
           DenseVector<Y> &y)
{
    assert(j>=basis.j0);
    assert(x.range()==basis.mra.rangeI(j+1));
    typedef typename X::ElementType T;

    basis.setLevel(j);
    DenseVector<Array<T> > sf, w;
    sf = transpose(basis.mra.M0) * x;
    w  = transpose(basis.M1) * x;
    concatenate(sf,w,y,y.firstIndex());
}

template <typename X, typename Y, Construction Cons>
void
reconstruct(const DenseVector<X> &x, 
            const Basis<typename X::ElementType,Primal,Interval,Cons> &basis, int j,
            DenseVector<Y> &y)
{
    typedef typename X::ElementType T;
    assert(j>=basis.j0);
    assert(x.range()==basis.mra.rangeI(j+1));
    basis.setLevel(j);
    y  = basis.mra.M0 * x(basis.mra.rangeI(j));
    y += basis.M1*x(_(basis.mra.rangeI(j).lastIndex()+1,basis.mra.rangeI(j+1).lastIndex()));
}

template <typename X, typename Y, Construction Cons>
void
reconstruct_(const DenseVector<X> &x,
             const Basis<typename X::ElementType,Dual,Interval,Cons> &basis_, int j,
             DenseVector<Y> &y)
{
    typedef typename X::ElementType T;
    assert(j>=basis_.j0);
    assert(x.range()==basis_.mra_.rangeI_(j+1));
    basis_.setLevel(j);
    y  = basis_.mra_.M0_ * x(basis_.mra_.rangeI_(j));
    y += basis_.M1_*x(_(basis_.mra_.rangeI_(j).lastIndex()+1,basis_.mra_.rangeI_(j+1).lastIndex()));
}

template <typename X, typename Y, Construction Cons>
void
fwt(const DenseVector<X> &x, 
    const Basis<typename X::ElementType,Dual,Interval,Cons> &basis_, int j,
    DenseVector<Y> &y)
{
    assert(j>=basis_.j0);
    assert(x.range()==basis_.mra_.rangeI_(j+1));
    
    y = x;
    for (int l=j; l>=basis_.j0; --l) {
        typename DenseVector<Y>::View yview = y(basis_.mra_.rangeI_(l+1));
        decompose(yview, basis_, l, yview);
    }
}

template <typename X, typename Y, Construction Cons>
void
ifwt(const DenseVector<X> &x, 
     const Basis<typename X::ElementType,Primal,Interval,Cons> &basis, int j,
     DenseVector<Y> &y)
{
    assert(j>=basis.j0);
    assert(x.range()==basis.mra.rangeI(j+1));
    
    y = x;
    for (int l=basis.j0; l<=j; ++l) {
        typename DenseVector<Y>::View yview = y(basis.mra.rangeI(l+1));
        DenseVector<Array<typename X::ElementType> > z = y(basis.mra.rangeI(l+1));
        reconstruct(z, basis, l, yview);
    }
}

} // namespace lawa

