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


#ifndef APPLICATIONS_FINANCE_SPACETIMECGMYOPERATOR1D_H
#define APPLICATIONS_FINANCE_SPACETIMECGMYOPERATOR1D_H 1


#include <lawa/integrals/integrals.h>
#include <lawa/operators/pdeoperators1d/pdeoperators1d.h>

namespace lawa {

template <typename T, typename Basis, typename CGMYOperator>
class SpaceTimeCGMYOperator1D{

public:

    typedef typename Basis::FirstBasisType Basis_t;
    typedef typename Basis::SecondBasisType Basis_x;

    typedef typename Basis_t::BSplineType PrimalSpline_t;
    typedef typename Basis_x::BSplineType PrimalSpline_x;
    typedef typename Basis_t::WaveletType PrimalWavelet_t;
    typedef typename Basis_x::WaveletType PrimalWavelet_x;

    typedef ConvectionOperator1D<T, Basis_t>   Convection_t;
    typedef IdentityOperator1D<T, Basis_t>        Reaction_t;
    typedef IdentityOperator1D<T, Basis_x>        Reaction_x;


    const Basis& basis;
    Convection_t d_t;
    Reaction_t   id_t;
    Reaction_x   id_x;
    const CGMYOperator &cgmy_x;


        PrimalSpline_t phi_t, d_phi_t;
        PrimalSpline_x phi_x;
        PrimalWavelet_t psi_t, d_psi_t;
        PrimalWavelet_x psi_x;

        Integral<T, Gauss, PrimalSpline_t, PrimalSpline_t>   integral_sfsf_t,
                                                           d_integral_sfsf_t;
        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x>   integral_sfsf_x;

        Integral<T, Gauss, PrimalSpline_t, PrimalWavelet_t>  integral_sfw_t,
                                                           d_integral_sfw_t;
        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>  integral_sfw_x;

        Integral<T, Gauss, PrimalWavelet_t, PrimalSpline_t>  integral_wsf_t,
                                                           d_integral_wsf_t;
        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>  integral_wsf_x;

        Integral<T, Gauss, PrimalWavelet_t, PrimalWavelet_t> integral_ww_t,
                                                           d_integral_ww_t;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x> integral_ww_x;


    public:
        SpaceTimeCGMYOperator1D(const Basis& _basis, const CGMYOperator &_cgmy_x);

        T getc() const;
        const Basis& getBasis() const;

        T                                                           // returns a(v,u)
        operator()(XType row_xtype_t, int j1_t, int k1_t,
                   XType row_xtype_x, int j1_x, int k1_x,
                   XType col_xtype_t,int j2_t, int k2_t,
                   XType col_xtype_x,int j2_x, int k2_x) const;
        T
        operator()(XType xtype_t, int j_t, int k_t,                 // returns a(u,u)
                   XType xtype_x, int j_x, int k_x) const;

        T
        operator()(const Index2D &row_index, const Index2D &col_index) const;


};

}    //namespace

#include <applications/finance/spacetimecgmyoperator1d.tcc>


#endif //APPLICATIONS_FINANCE_SPACETIMECGMYOPERATOR1D_H
