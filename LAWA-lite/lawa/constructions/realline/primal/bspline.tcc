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
#include <lawa/flensforlawa.h>

#include <lawa/math/math.h>
#include <extensions/extensions.h>

namespace lawa {

using namespace flens;

template <typename T>
    DenseVector<Array<T> >
    _bspline_mask(int d);

template <typename T>
BSpline<T,Primal,R,CDF>::BSpline(int _d)
    : d(_d), mu(d&1), l1(.5*(-d+mu)), l2(.5*(d+mu)),
      a(_bspline_mask<T>(d))
{
    assert(_d>0);
}

template <typename T>
BSpline<T,Primal,R,CDF>::BSpline(const MRA<T,Primal,R,CDF> &mra)
    : d(mra.d), mu(mra.d&1), l1(.5*(-d+mu)), l2(.5*(d+mu)), a(_bspline_mask<T>(d))
{
    assert(mra.d>0);
}

template <typename T>
BSpline<T,Primal,R,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,R,CDF>::operator()(T x, int j, long k, unsigned short deriv) const
{
    if (inner(x,support(j,k))) {
        T ret = T(0);
        x = pow2i<T>(j)*x-k - mu/2.;
        if (deriv==0) {
            x = fabs(x);
            for (int p=0; p<=ifloor(d/2.-x); ++p) {
                int sign = (p&1) ? -1 : 1;
                ret +=   sign * binomial(d, p) * pow(T(d/2.-x-p), d-1);
            }
            ret /= factorial(d - 1);
        } else {
            for (int p=0; p<=ifloor(d/2.-fabs(x)); ++p) {
                int sign = ( (p&1)==( (x>0)&&(deriv&1) ) ) ? 1 : -1;
                ret += sign * binomial(d, p)
                            * pow(T(d/2.-fabs(x)-p), d-1-deriv);
            }
            ret /= factorial(d-1-deriv);
        }
        // 2^(j*deriv) * 2^(j/2) * ret
        return pow2ih<T>(2*j*deriv+j)*ret;
    }
    return T(0);
}

template <typename T>
Support<T>
BSpline<T,Primal,R,CDF>::support(int j, long k) const
{
    return pow2i<T>(-j) * Support<T>(l1 + k, l2 + k);
}


template <typename T>
DenseVector<Array<T> >
BSpline<T,Primal,R,CDF>::singularSupport(int j, long k) const
{
    return linspace(support(j,k).l1, support(j,k).l2, d+1);
}

template <typename T>
T
BSpline<T,Primal,R,CDF>::tic(int j) const
{
    return pow2i<T>(-j);
}

template <typename T>
const DenseVector<Array<T> > &
BSpline<T,Primal,R,CDF>::mask() const
{
    return a;
}

//------------------------------------------------------------------------------

template <typename T>
BSpline<T,Primal,R,CDF>
N(int d)
{
    return BSpline<T,Primal,R,CDF>(d);
}

template <typename T>
DenseVector<Array<T> >
N1()
{
    DenseVector<Array<T> > ret(_(0, 1));
    ret = T(1),
            1;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N2()
{
    DenseVector<Array<T> > ret(_(-1, 1));
    ret = T(0.5),
            1,
            0.5;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N3()
{
    DenseVector<Array<T> > ret(_(-1, 2));
    ret = T(0.25),
            0.75,
            0.75,
            0.25;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N4()
{
    DenseVector<Array<T> > ret(_(-2, 2));
    ret = T(0.125),
            0.5,
            0.75,
            0.5,
            0.125;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N5()
{
    DenseVector<Array<T> > ret(_(-2, 3));
    ret = T(0.0625),
            0.3125,
            0.625,
            0.625,
            0.3125,
            0.0625;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N6()
{
    DenseVector<Array<T> > ret(_(-3, 3));
    ret = T(0.03125),
            0.1875,
            0.46875,
            0.625,
            0.46875,
            0.1875,
            0.03125;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N7()
{
    DenseVector<Array<T> > ret(_(-3, 4));
    ret = T(0.015625),
            0.109375,
            0.328125,
            0.546875,
            0.546875,
            0.328125,
            0.109375,
            0.015625;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N8()
{
    DenseVector<Array<T> > ret(_(-4, 4));
    ret = T(0.0078125),
            0.0625,
            0.21875,
            0.4375,
            0.546875,
            0.4375,
            0.21875,
            0.0625,
            0.0078125;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N9()
{
    DenseVector<Array<T> > ret(_(-4, 5));
    ret = T(0.00390625),
            0.03515625,
            0.140625,
            0.328125,
            0.4921875,
            0.4921875,
            0.328125,
            0.140625,
            0.03515625,
            0.00390625;
    return ret;
}

template <typename T>
DenseVector<Array<T> >
N10()
{
    DenseVector<Array<T> > ret(_(-5, 5));
    ret = T(0.001953125),
            0.01953125,
            0.087890625,
            0.234375,
            0.41015625,
            0.4921875,
            0.41015625,
            0.234375,
            0.087890625,
            0.01953125,
            0.001953125;
    return ret;
}

//------------------------------------------------------------------------------

template <typename T>
DenseVector<Array<T> >
_calculate_bspline_mask(int d)
{
    assert(d>=0);

    if (d==0) {
        DenseVector<Array<T> > res(1);
        res = T(1.);
        return res;
    }
    T factor = 1 << (d-1);
    int kappa = d & 1;

    int from = -(d+kappa)/2;
    int to   =  (d+kappa)/2;
    DenseVector<Array<T> > res(_(from,to));
    for (int i=from; i<=to; ++i) {
        res(i) = binomial(d,i-from) / factor;
    }
    return res;
}

template <typename T>
DenseVector<Array<T> >
_bspline_mask(int d)
{
    assert(d>=1);

    if (d>10) {
        return _calculate_bspline_mask<T>(d);
    }

    typedef DenseVector<Array<T> > (*NFunc)();
    static NFunc _NFunc[11] = {
        NULL,
        &N1,
        &N2,
        &N3,
        &N4,
        &N5,
        &N6,
        &N7,
        &N8,
        &N9,
        &N10
    };

    return _NFunc[d]();
}

} // namespace lawa

