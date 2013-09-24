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
#include <lawa/flensforlawa.h>

namespace flens {

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,Periodic,CDF>::RefinementMatrix(const BSpline<T,Side,Periodic,CDF> &spline)
      : band(Const<T>::R_SQRT2 * spline.mask())
{
}

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,Periodic,CDF>::RefinementMatrix(const Wavelet<T,Side,Periodic,CDF> &wavelet)
      : band(Const<T>::R_SQRT2 * wavelet.mask())
{
}


//------------------------------------------------------------------------------

template <typename X, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,Periodic,CDF> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y)
{
    typedef typename X::ElementType T;

    assert(alpha==1.);
    assert(x.engine().stride()==1);

    const DenseVector<Array<T> > &a = A.band;
    int lx = x.length();
    int la = a.length();
    
    if (transA==NoTrans) {
        const T *xp = x.engine().data();
        if (beta==0) {
            y.engine().resize(2*lx,x.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==2*lx);
        }
        int ly = y.length();
        if (la>ly) {
            // TODO: eliminate z (only needed since y can be a view and 
            // calculations are 0-based. Thus y can not be re-indexed.
            DenseVector<Array<T> > z = y;
            z.engine().changeIndexBase(0);
            for (int k=0; k<lx; ++k, ++xp) {
                int pos = (z.lastIndex() - (-a.firstIndex() % ly) + 1 + 2*k)%ly;
                for (int i=a.firstIndex(); i<=a.lastIndex(); ++i) {
                    if (pos>z.lastIndex()) {
                        pos = z.firstIndex();
                    }
                    z(pos++) += a(i) * *xp;
                }
            }
            y = z;
        } else {
            for (int k=0; k<lx; ++k, ++xp) {
                int mMin = a.firstIndex() + 2*k;
                int mMax = a.lastIndex() + 2*k;
                const T *abegin = a.engine().data(), *aend = abegin + a.length();
                T * ybegin = y.engine().data(), *yend = ybegin + y.length();
                if (mMin<0) {
                    cxxblas::axpy(-mMin  , *xp, abegin,      1, yend+mMin, 1);
                    cxxblas::axpy(la+mMin, *xp, abegin-mMin, 1, ybegin,    1);
                } else if (mMax>=ly) {
                    int p = mMax-ly+1;
                    cxxblas::axpy(p,    *xp, aend-p, 1, ybegin,     1);
                    cxxblas::axpy(la-p, *xp, abegin, 1, yend-la+p,   1);
                } else {
                    cxxblas::axpy(la,   *xp, abegin, 1, ybegin+mMin,1);
                }
            }
        }
    } else { // (transA==Trans)
        int la = a.length();
        if (beta==0) {
            y.engine().resize(lx/2,x.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==lx/2);
        }
        T *iter = y.engine().data();
        if (lx<=2*la-2) {
            const T *xbegin = x.engine().data(), *xend = xbegin + x.length();
            int shift = lx - (-a.firstIndex() % lx);
            for (int m=0; m<lx/2; ++m, ++iter) {
                const T *sigIter = xbegin+shift;
                *iter = 0;
                for (int i=a.firstIndex(); i<=a.lastIndex(); ++i, ++sigIter) {
                    if (sigIter==xend) {
                        sigIter = xbegin;
                    }
                    *iter += *sigIter * a(i);
                }
                if ((shift += 2) >= lx) {
                    shift = shift - lx;
                }
            }
        } else {
            DenseVector<Array<T> > ends(2*la-2);
            ends(_(1,la-1)) = x(_(x.firstIndex()+lx-la+1,x.firstIndex()+lx-1));
            ends(_(la,2*la-2)) = x(_(x.firstIndex(),x.firstIndex()+la-2));
            T *iter = y.engine().data();
            for (int m=0; m<lx/2; ++m, ++iter) {
                int kMin = a.firstIndex()+2*m;
                int kMax = a.lastIndex()+2*m;
                const T *abegin = a.engine().data();
                T dotValue;
                if ((kMin<0)) {
                    cxxblas::dot(la, abegin, 1, 
                                 ends.engine().data()+kMin+la-1, 1, dotValue);
                } else if ((kMax>=lx)) {
                    cxxblas::dot(la, abegin, 1, 
                                 ends.engine().data()+kMax-lx, 1, dotValue);
                } else {
                    cxxblas::dot(la, abegin, 1, 
                                 x.engine().data()+kMin, 1, dotValue);
                }
                *iter += dotValue;
            }
        }
    }
}

} // namespace flens

