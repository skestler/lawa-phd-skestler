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
#include <cmath>

#include <lawa/math/math.h>

namespace lawa {

template <typename T>
MRA<T,Primal,Interval,DKU>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), mu(d&1),
      l1((mu-d)/2), l2((mu+d)/2),
      l1_(l1-d_+1), l2_(l2+d_-1),
      l(l2_-(d_-d)), q(l+mu-1),
      min_j0(iceil<int>(log2(l2_+l2_-1)+1)),
      j0(std::max(j,min_j0)),
      _bc(2,0), _j(j0)
{
    assert(_j>=min_j0);
    assert((d+d_)%2==0);
    setLevel(_j);

    _alpha_();
    _beta_();

    _calcM0();
}

template <typename T>
MRA<T,Primal,Interval,DKU>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets (primal). ----------

template <typename T>
int
MRA<T,Primal,Interval,DKU>::cardI(int j) const
{
    assert(j>=min_j0);
    return d + pow2i<T>(j)-q-l + d - (_bc(0)+_bc(1));
}

template <typename T>
int
MRA<T,Primal,Interval,DKU>::cardIL(int /*j*/) const
{
    return d-_bc(0);
}

template <typename T>
int
MRA<T,Primal,Interval,DKU>::cardII(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j)-q-l;
}

template <typename T>
int
MRA<T,Primal,Interval,DKU>::cardIR(int /*j*/) const
{
    return d-_bc(1);
}

//--- ranges of whole, left, inner, right index sets (primal). -----------------

template <typename T>
Range<int>
MRA<T,Primal,Interval,DKU>::rangeI(int j) const
{
    assert(j>=min_j0);
    return Range<int>(l-d+_bc(0), pow2i<T>(j)-l+d-mu-_bc(1));
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,DKU>::rangeIL(int /*j*/) const
{
    return Range<int>(l-d+_bc(0), l-1);
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,DKU>::rangeII(int j) const
{
    assert(j>=min_j0);
    return Range<int>(l, pow2i<T>(j)-q-1);
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,DKU>::rangeIR(int j) const
{
    assert(j>=min_j0);
    return Range<int>(pow2i<T>(j)-q, pow2i<T>(j)-q+d-1-_bc(1));
}

template <typename T>
int
MRA<T,Primal,Interval,DKU>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Interval,DKU>::setLevel(int j)
{
    if (j!=_j) {
        assert(j>=min_j0);
        _j = j;
    }
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Primal,Interval,DKU>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;
}

//------------------------------------------------------------------------------

template <typename T>
void
MRA<T,Primal,Interval,DKU>::_alpha_()
{
    _Alpha_.engine().resize(_(1-l2, 2*l+l1-1), _(0,d-1));

    _Alpha_(_,0) = T(1.);      // (5.1.1)

    // 5.1.3
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    for (int r=1; r<d; ++r) {
        T tmp1 = 0.;
        for (int k=l1_; k<=l2_; ++k) {
            T tmp2 = 0.;
            for (int s=0; s<r; ++s) {
                tmp2 += binomial(r,s)*powii(k,r-s)*_Alpha_(0,s);
            }
            tmp1 += phi_.a_(k)*tmp2;
        }
        _Alpha_(0,r) = T(1.)/(pow2i<T>(r+1)-2) * tmp1;
    }

    // (5.1.2)
    for (int r=1; r<d; ++r) {
        // m<0
        for (int m=_Alpha_.firstRow(); m<0; ++m) {
            T tmp = 0.;
            for (int i=0; i<=r; ++i) {
                tmp += binomial(r,i)*powii(m,i)*_Alpha_(0,r-i);
            }
            _Alpha_(m,r) = tmp;
        }
        // m>0
        for (int m=1; m<=_Alpha_.lastRow(); ++m) {
            T tmp = 0.;
            for (int i=0; i<=r; ++i) {
                tmp += binomial(r,i)*powii(m,i)*_Alpha_(0,r-i);
            }
            _Alpha_(m,r) = tmp;
        }
    }
}

template <typename T>
void
MRA<T,Primal,Interval,DKU>::_beta_()
{
    if (d==1) {
        _Beta_.engine().resize(_(2*l+l1, 2*l+l1), _(0,0));
        return;
    }
    _Beta_.engine().resize(_(2*l+l1, 2*l+l2-2), _(0,d-1));

    // 3.2.31
    BSpline<T,Primal,R,CDF> phi(d);
    T factor = Const<T>::R_SQRT2;
    for (int r=0; r<d; ++r) {
        for (int m=_Beta_.firstRow(); m<=_Beta_.lastRow(); ++m) {
            T tmp = 0.;
            for (int q=iceil<int>((m-l2)/2.); q<l; ++q) {
                tmp += _Alpha_(q,r) * phi.a(m-2*q);
            }
            _Beta_(m,r) = factor * tmp;
        }
    }
}

template <typename T>
void
MRA<T,Primal,Interval,DKU>::_calcM0()
{
    GeMatrix<FullStorage<T,ColMajor> > Mj0(rangeI(min_j0+1), rangeI(min_j0)),
                                       Left(_(l-d,2*l+l2-2),_(l-d,l-1)), Right;

    for (int r=0; r<d; ++r) {
        Left(l-d+r, l-d+r) = Const<T>::R_SQRT2 * pow2i<T>(-r);
    }

    for (int m=l; m<=2*l+l1-1; ++m) {
        for (int k=0; k<d; ++k) {
            Left(m,l-d+k) = Const<T>::R_SQRT2 * pow2i<T>(-k) * _Alpha_(m,k);
        }
    }

    for (int m=2*l+l1; m<=2*l+l2-2; ++m) {
        for (int k=0; k<d; ++k) {
            Left(m,l-d+k) = _Beta_(m,k);
        }
    }
    Mj0(Left) = Left;
    arrow(Left,Right);
    Right.engine().changeIndexBase(Mj0.lastRow()-Right.numRows()+1,
                                   Mj0.lastCol()-Right.numCols()+1);
    Mj0(Right) = Right;
    
    Mj0.engine().changeIndexBase(1,1);

    BSpline<T,Primal,R,CDF> phi(d);
    for (int c=l, r=2*l+l1; c<=pow2i<T>(min_j0)-q-1; ++c, r+=2) {
        Mj0(_(r,r+phi.a.length()-1),c) = Const<T>::R_SQRT2 * phi.a; 
    }

    for (int c=l, r=2*l+l1; c<=pow2i<T>(min_j0)-q-1; ++c, r+=2) {
        Mj0(_(r,r+phi.a.length()-1),c) = Const<T>::R_SQRT2 * phi.a; 
    }
    
    M0 = RefinementMatrix<T,Interval,DKU>(d,d,Mj0,min_j0);
}

} // namespace lawa

