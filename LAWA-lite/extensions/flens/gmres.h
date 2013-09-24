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

#ifndef EXTENSIONS_FLENS_GMRES_H
#define EXTENSIONS_FLENS_GMRES_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

using namespace flens;

template <typename A>
struct _gmres
{
    typedef typename A::ElementType T;
};

template <typename MA, typename VX, typename VB>
int
gmres(const MA &A, VX &x, const VB &b,
      typename _gmres<MA>::T tol = std::numeric_limits<typename _gmres<MA>::T>::epsilon(),
      long maxIterations = -1);

template <typename MA, typename VX, typename VB>
int
gmresm(const MA &A, VX &x, const VB &b,
       typename _gmres<MA>::T tol=std::numeric_limits<typename _gmres<MA>::T>::epsilon(),
       int restart=20, long maxIterations=-1);

template <typename Prec, typename MA, typename VX, typename VB>
int
pgmres(const Prec &Pr, const MA &A, VX &x, const VB &b,
       typename _gmres<MA>::T tol = std::numeric_limits<typename _gmres<MA>::T>::epsilon(),
       long maxIterations = -1);

template <typename Prec, typename MA, typename VX, typename VB>
int
pgmresm(const Prec &Pr, const MA &A, VX &x, const VB &b,
        typename _gmres<MA>::T tol=std::numeric_limits<typename _gmres<MA>::T>::epsilon(),
        int restart=20, long maxIterations=-1);

//-- trait specialization for FLENS types ------------------------------------

template <typename I>
struct _gmres<Matrix<I> >
{
    typedef typename Matrix<I>::ElementType T;
};

} // namespace lawa

#include <extensions/flens/gmres.tcc>

#endif // EXTENSIONS_FLENS_GMRES_H

