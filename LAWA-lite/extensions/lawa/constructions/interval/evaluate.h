/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_EVALUATE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_EVALUATE_H 1

#include <lawa/constructions/basis.h>
#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <FunctionSide Side, Construction Cons, typename X>
    typename X::ElementType
    evaluate(const MRA<typename X::ElementType,Side,Interval,Cons> &mra,
             int j, const DenseVector<X> &coeffs, typename X::ElementType x, 
             int deriv);

template <FunctionSide Side, Construction Cons, typename X>
    typename X::ElementType
    evaluate(const Basis<typename X::ElementType,Side,Interval,Cons> &basis,
             int J, const DenseVector<X> &coeffs, typename X::ElementType x, 
             int deriv);

} // namespace lawa

#include <lawa/constructions/interval/evaluate.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_EVALUATE_H

