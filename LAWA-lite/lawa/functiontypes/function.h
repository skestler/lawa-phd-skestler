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

#ifndef LAWA_FUNCTIONTYPES_FUNCTION_H
#define LAWA_FUNCTIONTYPES_FUNCTION_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

using namespace flens;

template <typename T>
struct Function
{
    Function(T (*_f)(T));

    Function(T (*_f)(T), const DenseVector<Array<T> > &_singularPoints);

    T
    operator()(T x) const;

    T (*f)(T);
    const DenseVector<Array<T> > &singularPoints;
};

} // namespace lawa

#include <lawa/functiontypes/function.tcc>

#endif //LAWA_FUNCTIONTYPES_FUNCTION_H

