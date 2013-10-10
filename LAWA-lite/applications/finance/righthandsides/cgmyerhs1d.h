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

#ifndef APPLICATIONS_FINANCE_RIGHTHANDSIDES_CGMYERHS1D_H
#define APPLICATIONS_FINANCE_RIGHTHANDSIDES_CGMYERHS1D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/operators/deltas.h>
#include <lawa/settings/settings.h>
#include <applications/finance/options/options.h>
#include <applications/finance/righthandsides/optionrhs1d.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <typename T, typename Basis1D>
struct OptionRHS1D<T, Put, CGMYe, Basis1D>
{
    OptionRHS1D(const OptionParameters1D<T,Put> &_optionparameters,
                const ProcessParameters1D<T,CGMYe> &_processparameters,
                const Basis1D &_basis, T _R1=0., T _R2=1., bool _excessToPayoff=true);

    T
    operator()(XType xtype, int j, int k) const;

    T
    operator()(T time, XType xtype, int j, int k) const;

    T
    operator()(const Index1D &lambda) const;

    T
    operator()(T time, const Index1D &lambda) const;

    const OptionParameters1D<T,Put>     &optionparameters;
    const ProcessParameters1D<T,CGMYe>   &processparameters;
    const Basis1D                       &basis;
    Kernel<T,CGMY>                      kernel;
    T                                   R1, R2;
    T                                   OneDivSqrtR2pR1, OneDivR2pR1, R1DivR1pR2;
    bool                                excessToPayoff;
};

}   // namespoce lawa

#include <applications/finance/righthandsides/cgmyerhs1d.tcc>

#endif  // APPLICATIONS_FINANCE_RIGHTHANDSIDES_CGMYERHS1D_H
