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

#ifndef EXTENSIONS_FLENS_POWERMETHOD_H
#define EXTENSIONS_FLENS_POWERMETHOD_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

template <typename A>
struct _powerMethod
{
    typedef typename A::ElementType T;
};

template <typename MatA, typename VecX>
    void
    powerMethod(const MatA &A, typename MatA::ElementType tol,
                typename MatA::ElementType &lambda, VecX &y);

template <typename MatA, typename VecX>
    void
    inversePowerMethod(const MatA &A, typename MatA::ElementType tol,
                       typename MatA::ElementType &lambda, VecX &y);

/*

template <typename Prec, typename MatA, typename VecX>
    void
    powerMethod(const Prec &P, const MatA &A, typename MatA::ElementType tol,
                typename MatA::ElementType &lambda, VecX &y);

//--- SymmetricMatrix ---

template <typename I, typename VecX>
    void
    inversePowerMethod(const SymmetricMatrix<I> &A,
                       typename _powerMethod<I>::T tol,
                       typename _powerMethod<I>::T &lambda, VecX &y);

template <typename Prec, typename I, typename VecX>
    void
    inversePowerMethod(const Prec &P, const SymmetricMatrix<I> &A,
                       typename _powerMethod<I>::T tol,
                       typename _powerMethod<I>::T &lambda, VecX &y);

//--- Matrix ---

template <typename I, typename VecX>
    void
    inversePowerMethod(const Matrix<I> &A, typename _powerMethod<I>::T tol,
                       typename _powerMethod<I>::T &lambda, VecX &y);

template <typename Prec, typename I, typename VecX>
    void
    inversePowerMethod(const Prec &P, const Matrix<I> &A, typename _powerMethod<I>::T tol,
                       typename _powerMethod<I>::T &lambda, VecX &y);

template <typename MatA>
    typename MatA::ElementType
    condition(const MatA &A, typename MatA::ElementType tol=1e-12);


template <typename Prec, typename MatA>
    typename MatA::ElementType
    condition(const Prec &P, const MatA &A, typename MatA::ElementType tol=1e-12);

*/
//-- trait specialization for FLENS types --------------------------------------

template <typename I>
struct _powerMethod<Matrix<I> >
{
    typedef typename Matrix<I>::ElementType T;
};

} // namespace lawa

#include <extensions/flens/powermethod.tcc>

#endif // EXTENSIONS_FLENS_POWERMETHOD_H

