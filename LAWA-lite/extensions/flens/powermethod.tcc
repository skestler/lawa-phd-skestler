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

#include <cstdlib>

#include <extensions/flens/cg.h>
#include <extensions/flens/gmres.h>
#include <extensions/flens/vectornorm.h>

namespace lawa {

template <typename MatA, typename VecX>
void
powerMethod(const MatA &A, typename MatA::ElementType tol,
            typename MatA::ElementType &lambda, VecX &y)
{
    typedef typename MatA::ElementType T;
    int maxIterations = std::numeric_limits<int>::max();

    lambda = 0.;
    T lambdaOld = 0.;
    VecX z = y;
    for (int i=1; i<=maxIterations; ++i) {
        lambdaOld = lambda;
        y = A*z;
        lambda = norm<lInfinity>(y);
        if (fabs(lambda-lambdaOld)<tol*lambda) {
            lambda = (y*z)/(z*z);
            y = z;
            return;
        }
        z = ((T)1./lambda) * y;
    }
}

template <typename MatA, typename VecX>
void
inversePowerMethod(const MatA &A, typename MatA::ElementType tol,
                   typename MatA::ElementType &lambda, VecX &y)
{
    typedef typename MatA::ElementType T;
    int maxIterations = std::numeric_limits<int>::max();

    VecX v, q;
    for (int i=1; i<=maxIterations; ++i) {
        v = (1./norm<l2>(y)) * y;
        gmresm(A,y,v,tol);
        lambda = v*y;
        if (norm<l2>(q = -lambda*v+y) <= tol*fabs(lambda)) {
            lambda = 1./lambda;
            return;
        }
    }
}

/*

template <typename Prec, typename MatA, typename VecX>
void
powerMethod(const Prec &P, const MatA &A, typename MatA::ElementType tol,
            typename MatA::ElementType &lambda, VecX &y)
{
    typedef typename MatA::ElementType T;
    int maxIterations = std::numeric_limits<int>::max();

    lambda = 0.;
    T lambdaOld = 0.;
    VecX z = y, w;
    for (int i=1; i<=maxIterations; ++i) {
        lambdaOld = lambda;
        y = P*z;
        w = A*y;
        y = transpose(P)*w;
        lambda = norm<lInfinity>(y);
        if (absolute(lambda-lambdaOld)<tol*lambda) {
            lambda = (y*z)/(z*z);
            y = z;
            return;
        }
        z = y/lambda;
    }
}

//--- SymmetricMatrix ---

template <typename I, typename VecX>
void
inversePowerMethod(const SymmetricMatrix<I> &A, typename _powerMethod<I>::T tol,
                   typename _powerMethod<I>::T &lambda, VecX &y)
{
    int maxIterations = std::numeric_limits<int>::max();

    VecX v, q;
    for (int i=1; i<=maxIterations; ++i) {
        v = y / norm<l2>(y);
        lawa::cg(A,y,v,tol);
        lambda = v*y;
        if (norm<l2>(q = -lambda*v+y) <= tol*absolute(lambda)) {
            lambda = 1./lambda;
            return;
        }
    }
}

template <typename Prec, typename I, typename VecX>
void
inversePowerMethod(const Prec &P, const SymmetricMatrix<I> &A, typename _powerMethod<I>::T tol,
                   typename _powerMethod<I>::T &lambda, VecX &y)
{
    int maxIterations = std::numeric_limits<int>::max();

    VecX v, q = y, w = y;
    for (int i=1; i<=maxIterations; ++i) {
        v = y / norm<l2>(y);
        lawa::cg(P,w,v,tol);
        lawa::pcg(P,A,q,w,tol);
        lawa::cg(P,y,q,tol);
        lambda = v*y;
        if (norm<l2>(q = -lambda*v+y) <= tol*absolute(lambda)) {
            lambda = 1./lambda;
            return;
        }
    }
}

//--- Matrix ---

template <typename I, typename VecX>
void
inversePowerMethod(const Matrix<I> &A, typename _powerMethod<I>::T tol,
                   typename _powerMethod<I>::T &lambda, VecX &y)
{
    int maxIterations = std::numeric_limits<int>::max();
    VecX v, q;
    for (int i=1; i<=maxIterations; ++i) {
        v = (1./norm<l2>(y)) * y;
        gmres(A,y,v,tol);
        lambda = v*y;
        if (norm<l2>(q = -lambda*v+y) <= tol*absolute(lambda)) {
            lambda = 1./lambda;
            return;
        }
    }
}

template <typename Prec, typename I, typename VecX>
void
inversePowerMethod(const Prec &P, const Matrix<I> &A, typename _powerMethod<I>::T tol,
                   typename _powerMethod<I>::T &lambda, VecX &y)
{
    int maxIterations = std::numeric_limits<int>::max();

    VecX v, q = y, w = y;
    for (int i=1; i<=maxIterations; ++i) {
        v = y / norm<l2>(y);
        lawa::gmres(transpose(P),w,v,tol);
        lawa::pgmres(P,A,q,w,tol);
        lawa::gmres(P,y,q,tol);
        lambda = v*y;
        if (norm<l2>(q = -lambda*v+y) <= tol*absolute(lambda)) {
            lambda = 1./lambda;
            return;
        }
    }
}

template <typename MatA>
typename MatA::ElementType
condition(const MatA &A, typename MatA::ElementType tol)
{
    typedef typename MatA::ElementType T;
    DenseVector<Array<T> > ev(dim(A));
    ev = T(1);
    T lambdaMax, lambdaMin;
    powerMethod(A, tol, lambdaMax, ev);
    ev = T(1);
    inversePowerMethod(A, tol, lambdaMin, ev);

    assert(lambdaMin!=0);
    return absolute(lambdaMax) / absolute(lambdaMin);
}

template <typename Prec, typename MatA>
typename MatA::ElementType
condition(const Prec &P, const MatA &A, typename MatA::ElementType tol)
{
    typedef typename MatA::ElementType T;
    DenseVector<Array<T> > ev(dim(A));
    ev = T(1);
    T lambdaMax, lambdaMin;
    powerMethod(P, A, tol, lambdaMax, ev);
    ev = T(1);
    inversePowerMethod(P, A, tol, lambdaMin, ev);

    assert(lambdaMin!=0);
    return absolute(lambdaMax) / absolute(lambdaMin);
}

*/

} // namespace lawa

