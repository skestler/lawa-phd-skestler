/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef APPLICATIONS_FINANCE_CGMYOPERATOR1D_FAST_H
#define APPLICATIONS_FINANCE_CGMYOPERATOR1D_FAST_H 1

#include <applications/finance/cgmyutils.h>

namespace lawa {


template <typename T, typename Basis>
class CGMYOperator1D_Fast{

    public:

        const Basis& basis;
        const T diffusion, convection, reaction;
        const T C,G,M,Y;

        CGMYUtils<T> cgmy;

        typedef typename Basis::BSplineType PrimalSpline;
        typedef typename Basis::WaveletType PrimalWavelet;

        PrimalSpline phi, d_phi, delta_phi;
        PrimalWavelet psi, d_psi, delta_psi;

        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, d_integral_sfsf, dd_integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet> integral_sfw, d_integral_sfw,  dd_integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline> integral_wsf, d_integral_wsf,  dd_integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, d_integral_ww,   dd_integral_ww;

    public:
        CGMYOperator1D_Fast(const Basis& _basis, T _diffusion, T _convection, T _reaction, T _C, T _G, T _M, T _Y);

        T getc() const;
        const Basis& getBasis() const;

        T
        operator()(XType xtype1, int j1, int k1,
                   XType xtype2, int j2, int k2) const;

        T
        operator()(XType xtype1, int j1, int k1,
                   XType xtype2, int j2, int k2, T R1, T R2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;

        GeMatrix<FullStorage<T,ColMajor> >
        computeDeltas(XType xtype, int j, int k) const;


};

}    //namespace lawa

#include <applications/finance/cgmyoperator1d_fast.tcc>


#endif // APPLICATIONS_FINANCE_CGMYOPERATOR1D_FAST_H
