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

#ifndef APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_RPLUS1D_H
#define APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_RPLUS1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

/*
 * Reference solutions u and corresponding righthand sides for second order PDEs
 * with constant coefficients on R:
 *       - diffusion * u'' + convection * u' + reaction * u = f
 */

template<typename T>
struct RefSols_PDE_RPlus1D
{
    static int nr;

    static T reaction, convection, diffusion;

    static DenseVector<Array<T> > sing_pts;

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas;
    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > H1_deltas;

    static void
    setExample(int _nr, T _reaction, T _convection, T diffusion);

    static T
    exact(T x, int deriv);

    static T
    u(T x);

    static T
    d_u(T x);

    static T
    rhs(T x);

    // required for post-processing
    static T
    H1_rhs(T x);

    static T
    H1norm();

};


} // namespace lawa

#include <applications/unbounded_domains/referencesolutions/refsols_pde_rplus1d.tcc>


#endif // APPLICATIONS_UNBOUNDED_DOMAINS_REFERENCESOLUTIONS_REFSOLS_PDE_RPLUS1D_H
