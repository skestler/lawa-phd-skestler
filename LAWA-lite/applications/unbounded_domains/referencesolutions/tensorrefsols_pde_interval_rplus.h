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

#ifndef APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_INTERVAL_RPLUS_H
#define APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_INTERVAL_RPLUS_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template<typename T>
struct TensorRefSols_PDE_Interval_RPlus
{
    static int nr;
    static T reaction;
    static T convection_x;
    static T convection_y;
    static T diffusion_y;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y;   //aligned singularities
    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas_x, deltas_y;
    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > H1_deltas_x, H1_deltas_y;

    // diffusion_x is assumed to 1!!
    static void
    setExample(int _nr, T _reaction, T _convection_x, T _convection_y, T _diffusion_y);

    static T
    exact(T x, T y);

    static T
    exact_dx(T x, T y);

    static T
    exact_dy(T x, T y);

    static T
    exact_x(T x);

    static T
    exact_x(T x, int deriv_x);

    static T
    exact_y(T y);

    static T
    exact_y(T y, int deriv_y);

    static T
    rhs_x(T x);

    static T
    H1_rhs_x(T x);

    static T
    rhs_y(T y);

    static T
    H1_rhs_y(T y);

    static T
    H1norm();

    static T
    Energynorm();
};

}   // namespace lawa

#include <applications/unbounded_domains/referencesolutions/tensorrefsols_pde_interval_rplus.tcc>

#endif  // APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_TENSORREFSOLS_PDE_INTERVAL_RPLUS_H
