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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_LOCALWEIGHTING2D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_LOCALWEIGHTING2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T, typename Basis>
struct LocalWeighting2D
{
    static T            left_x1;
    static T            right_x1;
    static T            left_x2;
    static T            right_x2;
    static T            RightmLeft_x1;
    static T            RightmLeft_x2;

    static const Basis  *basis;

    static int          weight_type;

    static void
    setDomain(T _left_x1, T _right_x1, T _left_x2, T _right_x2);

    static void
    setBasis(const Basis *_basis);

    static void
    setWeightType(int _weight_type);

    static T
    weight(const Index2D &index);
};

}   // namespace lawa

#include <applications/finance/initialconditions/localweighting2d.tcc>

#endif // APPLICATIONS_FINANCE_INITIALCONDITIONS_LOCALWEIGHTING2D_H
