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

#ifndef EXTENSIONS_FLENS_VECTORNORM_H
#define EXTENSIONS_FLENS_VECTORNORM_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

enum NormType { l1, l2, lInfinity};

template <NormType N, typename T>
struct NormTraits
{
    typedef T Type;
};

template <NormType N, typename X>
    typename NormTraits<N,typename X::ElementType>::Type
    norm(const DenseVector<X> &x);

template <NormType N, typename T>
    typename NormTraits<N,T>::Type
    norm(const Array<T> &x);

} // namespace flens

#include <extensions/flens/vectornorm.tcc>

#endif // FEXTENSIONS_FLENS_VECTORNORM_H
