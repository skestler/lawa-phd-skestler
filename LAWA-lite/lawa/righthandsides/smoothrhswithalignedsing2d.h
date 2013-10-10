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

#ifndef LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H
#define LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/functiontypes/function2d.h>
#include <lawa/integrals/integral2d.h>

namespace lawa {

template<typename T, typename Basis2D, QuadratureType Quad>
struct SmoothRHSWithAlignedSing2D
{
    SmoothRHSWithAlignedSing2D(const Basis2D& _basis, const Function2D<T>& _F,
                               int order, unsigned short _derivx=0, unsigned short _derivy=0);

    T
    operator()(XType xtype_x, int j_x, long k_x,
               XType xtype_y, int j_y, long k_y) const;

    T
    operator()(const Index2D &index) const;

    const Basis2D &basis;
    Integral2D<Quad, typename Basis2D::FirstBasisType,
                     typename Basis2D::SecondBasisType> integral2d;
    const unsigned short derivx, derivy;


};


}    //namespace lawa

#include <lawa/righthandsides/smoothrhswithalignedsing2d.tcc>


#endif //LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H 1

