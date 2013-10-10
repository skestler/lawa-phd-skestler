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
 
#ifndef APPLICATIONS_FINANCE_KERNELS_CGMYKERNEL_H
#define APPLICATIONS_FINANCE_KERNELS_CGMYKERNEL_H 1

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <applications/finance/kernels/kernel.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <typename T>
struct Kernel<T,CGMY>
{
    typedef boost::math::policies::policy<boost::math::policies::digits2<64> > my_prec_policy;

    Kernel(const ProcessParameters1D<T,CGMY>  &_params);

    Kernel(const ProcessParameters1D<T,CGMYe> &_params);

    ProcessParameters1D<T,CGMY> params;
    const T C, G, M, Y;
    T c3,c4,c5,c6;
    T ExpXmOnemX_k_pos, ExpXmOnemX_k_neg, ExpXmOnemX_k;
    T ExpXmOne_k1_pos, ExpXmOne_k1_neg;
    T ExpXmOne_k2_pos, ExpXmOne_k2_neg;
    T constants[10];

    T powM_Ym5, powM_Ym4, powM_Ym3, powM_Ym2, powM_Ym1, powM_Y;
    T powG_Ym5, powG_Ym4, powG_Ym3, powG_Ym2, powG_Ym1, powG_Y;
    T CdivY, Cdiv1mY, Mdiv1mY, Gdiv1mY, OnedivSix;

    my_prec_policy Pol;

    T TailIntegralMoments(T x, int l) const;

    T nthTailIntegral(T x, int n, AntiDerivativeType type) const;


    T
    operator()(T x) const { return this->SecondTailIntegral(-x); }

    T
    FirstTailIntegral(T x) const;

    T
    FirstTailFirstExpMomIntegral(T x) const;

    T
    SecondTailIntegral(T x) const;

    T
    SecondTailFirstExpMomIntegral(T x) const;

    T
    ThirdTailIntegral(T x) const;

    T
    ThirdTailFirstExpMomIntegral(T x) const;

    T
    ForthTailIntegral(T x) const;

    T
    ForthTailFirstExpMomIntegral(T x) const;

    T
    ForthTailIntegral_test(T x) const;

    T
    FifthTailIntegral(T x) const;

    T
    SixthTailIntegral(T x) const;

};

}   // namespace lawa

#include <applications/finance/kernels/cgmykernel.tcc>

#endif  // APPLICATIONS_FINANCE_KERNELS_CGMYKERNEL_H

