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

#ifndef APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_REALLINE3D_H
#define APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_REALLINE3D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template<typename T>
struct TensorRefSols_PDE_Realline3D
{
    static int nr;
    static T c;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y, sing_pts_z;   //aligned singularities

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas_x, deltas_y, deltas_z;

    static void
    setExample(int _nr, T _c);

    static T
    exact(T x, T y, T z);

    static T
    exact_dx(T x, T y, T z);

    static T
    exact_dy(T x, T y, T z);

    static T
    exact_dz(T x, T y, T z);

    static T
    exact_x(T x);

    static T
    exact_x(T x, int deriv_x);

    static T
    exact_y(T y);

    static T
    exact_y(T y, int deriv_y);

    static T
    exact_z(T z);

    static T
    exact_z(T z, int deriv_z);

    static T
    rhs_x(T x);

    static T
    rhs_y(T y);

    static T
    rhs_z(T z);

    static T
    H1norm();
};

}   // namespace lawa

#include <applications/unbounded_domains/referencesolutions/tensorrefsols_pde_realline3d.tcc>

#endif  // APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_REALLINE3D_H
