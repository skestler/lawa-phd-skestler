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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_DIAGONALPRECONDITIONERADAPTIVEOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_DIAGONALPRECONDITIONERADAPTIVEOPERATOR_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator>
struct DiagonalPreconditionerAdaptiveOperator
{


    DiagonalPreconditionerAdaptiveOperator(AdaptiveOperator &_A);

    T
    operator()(const Index &index);

    AdaptiveOperator &A;

};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/diagonalpreconditioneradaptiveoperator.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_DIAGONALPRECONDITIONERADAPTIVEOPERATOR_H

