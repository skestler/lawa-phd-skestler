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

#ifndef LAWA_SPY_H
#define LAWA_SPY_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

using namespace flens;

template <typename T>
    void
    spy(const SparseGeMatrix<CRS<T,CRS_General> > &A, const char* filename,
        bool absoluteValues = true, T eps = std::numeric_limits<T>::epsilon());

template <typename I>
    void
    spy (const Matrix<I> &A, const char* filename, bool absoluteValues = true,
         typename I::ElementType eps = std::numeric_limits<typename I::ElementType>::epsilon());

} // namespace lawa

#include <extensions/flens/spy.tcc>

#endif // LAWA_SPY_H

