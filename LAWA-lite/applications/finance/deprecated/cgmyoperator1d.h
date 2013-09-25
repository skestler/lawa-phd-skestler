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

#ifndef APPLICATIONS_FINANCE_CGMYOPERATOR1D_H
#define APPLICATIONS_FINANCE_CGMYOPERATOR1D_H 1

#include <applications/finance/cgmyutils.h>

namespace lawa {

template <typename T, typename Basis>
class CGMYOperator1D
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

    public:
        CGMYOperator1D(const Basis &basis, T _diffusion, T _convection, T _reaction,
                              T C, T G, T M, T Y, int _order=20, int _n=10, T _sigma=0.3);

        T
        getc() const;

        const Basis&
        getBasis() const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;

        T
        operator()(XType xtype_row, int j_row, int k_row, XType xtype_col, int j_col, int k_col, const CGMYUtils<T> &cgmy_comp) const;

    private:
        const Basis &basis;
        T diffusion, convection, reaction;

        PrimalSpline phi, d_phi;
        PrimalWavelet psi, d_psi;
        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, d_integral_sfsf, dd_integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet> integral_sfw, d_integral_sfw,  dd_integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline> integral_wsf, d_integral_wsf,  dd_integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, d_integral_ww,   dd_integral_ww;

        CGMYUtils<T> cgmy, cgmy_adjoint;

        int n;                //parameters for composite variable order Gauss-Legendre quadrature
        T sigma;
        T omega;
        T mu;
        int order;            //fixed Legendre order for outer integral

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;


        T
        _integrate_dpsi_vs_int_dpsi_k(T a, T b, XType xtype_row, int j_row, int k_row,
                                                   XType xtype_col, int j_col, int k_col,
                                                   const CGMYUtils<T> &cgmy_comp) const;

        T
        _nonsingular_quadrature_dpsi_vs_CGMYkernel(T x_ast, T a, T b,
                                                   XType xtype_col, int j_col, int k_col,
                                                   const CGMYUtils<T> &cgmy_comp) const;

        T
        _singular_quadrature_dpsi_vs_CGMYkernel(T x_ast, T a, T b,
                                                XType xtype_col, int j_col, int k_col,
                                                bool cutoff,
                                                const CGMYUtils<T> &cgmy_comp) const;

};

}    //namespace lawa

#include <applications/finance/cgmyoperator1d.tcc>

#endif    //LAWA_OPERATORS_INTEGRALOPERATORS_CGMYOPERATOR1D_H
