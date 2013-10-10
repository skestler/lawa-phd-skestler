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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_H 1

namespace lawa {

template <typename T>
    T
    _constant_bspline_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _constant_wavelet_inner_evaluator0(T x, unsigned short deriv);

}   // namespace lawa

#include <lawa/constructions/interval/multi/_constant_evaluator.tcc>

#endif //LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_H
