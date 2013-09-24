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

#ifndef EXTENSIONS_GMPMATH_POW2_H
#define EXTENSIONS_GMPMATH_POW2_H 1

#include <gmpxx.h>
#include <lawa/aux/aux.h>

namespace lawa {

template <typename T>
    typename GMPReal<T>::Type
    pow2i(int expo);

template <typename T>
    typename GMPReal<T>::Type
    pow2ih(int expo);

} // namespace lawa

#include <extensions/gmpmath/pow2.tcc>

#endif //  EXTENSIONS_GMPMATH_POW2_H

