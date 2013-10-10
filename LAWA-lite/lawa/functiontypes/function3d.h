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

#ifndef LAWA_FUNCTIONTYPES_FUNCTION3D_H
#define LAWA_FUNCTIONTYPES_FUNCTION3D_H 1

namespace lawa {

using namespace flens;

template<typename T>
struct Function3D
{
public:
    Function3D(T (*_f)(T,T,T), const DenseVector<Array<T> > &_singularPts_x,
                             const DenseVector<Array<T> > &_singularPts_y,
                             const DenseVector<Array<T> > &_singularPts_z);

    T
    operator()(T x, T y, T z) const;

    const DenseVector<Array<T> > singularPts_x; //x-aligned singularities
    const DenseVector<Array<T> > singularPts_y; //y-aligned singularities
    const DenseVector<Array<T> > singularPts_z; //y-aligned singularities
    T (*f)(T,T,T);
};

} // namespace lawa

#include <lawa/functiontypes/function3d.tcc>

#endif // LAWA_FUNCTIONTYPES_FUNCTION3D_H

