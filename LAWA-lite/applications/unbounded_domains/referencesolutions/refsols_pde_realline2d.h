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

#ifndef APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_REALLINE2D_H
#define APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_REALLINE2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template<typename T>
struct RefSols_PDE_Realline2D
{
    static int nr;
    static T c;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y;

    static void
    setExample(int _nr, T _c);

    static T
    exact(T x, T y);

    static T
    minus_exact(T x, T y);

    static T
    exact_dx(T x, T y);

    static T
    exact_dy(T x, T y);

    static T
    rhs(T x, T y);

    static T
    exact(T x, T y, int deriv_x, int deriv_y);

    static T
    H1norm();

};

}   // namespace lawa

#include <applications/unbounded_domains/referencesolutions/refsols_pde_realline2d.tcc>

#endif  // APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_REALLINE2D_H
