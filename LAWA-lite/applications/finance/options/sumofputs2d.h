/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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


#ifndef APPLICATIONS_FINANCE_OPTIONS_SUMOFPUTS2D_H
#define APPLICATIONS_FINANCE_OPTIONS_SUMOFPUTS2D_H 1

#include <boost/math/distributions/normal.hpp>

#include <applications/finance/options/option.h>
#include <applications/finance/options/optionparameters2d.h>
#include <applications/finance/processes/processes.h>
#include <applications/finance/fourierpricing/fourierpricing.h>

namespace lawa {

template <typename T>
struct Option2D<T,SumOfPuts>
{
    Option2D(void);

    Option2D(const OptionParameters2D<T,SumOfPuts> &_optionparameters);

    T
    payoff(T S1, T S2) const;

    T
    payoff_log(T x1, T x2) const;

    T
    value(const ProcessParameters2D<T,BlackScholes2D> &processparameters, T S1, T S2, T t) const;

    T
    value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters, T S1, T S2, T t);

    OptionParameters2D<T,SumOfPuts> optionparameters;

    OptionParameters1D<T,Put>       optionparameters1, optionparameters2;
    Option1D<T,Put>                 option1, option2;
};

} // namespace lawa

#include <applications/finance/options/sumofputs2d.tcc>

#endif // APPLICATIONS_FINANCE_OPTIONS_SUMOFPUTS2D_H
