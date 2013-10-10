/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Primal,R,CDF>
    : public BasisFunction<_T,Primal,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    BSpline(int _d);

    BSpline(const MRA<T,Primal,R,CDF> &mra);

    virtual
    ~BSpline();

    T
    operator()(T x, int j, long k, unsigned short deriv) const;

    Support<T>
    support(int j, long k) const;

    DenseVector<Array<T> >
    singularSupport(int j, long k) const;

    T
    tic(int j) const;

    const DenseVector<Array<T> > &
    mask() const;

    const int d, mu;
    const int l1, l2;
    const DenseVector<Array<T> > a;
};

//------------------------------------------------------------------------------

template <typename T>
BSpline<T,Primal,R,CDF>
N(int d);

template <typename T>
DenseVector<Array<T> >
N1();

template <typename T>
DenseVector<Array<T> >
N2();

template <typename T>
DenseVector<Array<T> >
N3();

template <typename T>
DenseVector<Array<T> >
N4();

template <typename T>
DenseVector<Array<T> >
N5();

template <typename T>
DenseVector<Array<T> >
N6();

template <typename T>
DenseVector<Array<T> >
N7();

template <typename T>
DenseVector<Array<T> >
N8();

template <typename T>
DenseVector<Array<T> >
N9();

template <typename T>
DenseVector<Array<T> >
N10();

} // namespace lawa

#include <lawa/constructions/realline/primal/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H

