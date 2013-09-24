/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY2D_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

//template <typename T, typename Index, typename Basis, typename MA, typename BilinearForm>
//Coefficients<Lexicographical,T,Index>
//APPLY(MA &A, const Coefficients<Lexicographical,T,Index> &v, int k);

template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
APPLY_Helmholtz2D(MA &A, const Coefficients<Lexicographical,T,Index2D> &v, int k);


}    //namespace lawa

#include <lawa/methods/adaptive/algorithms/apply2d.tcc>

#endif //  LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY2D_H

