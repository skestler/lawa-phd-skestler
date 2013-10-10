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

#ifndef LAWA_OPERATORS_DELTAS_H
#define LAWA_OPERATORS_DELTAS_H 1

#include <lawa/settings/enum.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/constructions/support.h>


namespace lawa {

template <typename T, typename Basis>
GeMatrix<FullStorage<T,ColMajor> >
computeDeltas(const Basis &basis, int j, long k, XType e);

}    //namespace lawa

#include <lawa/operators/deltas.tcc>

#endif    //LAWA_OPERATORS_DELTAS_H
