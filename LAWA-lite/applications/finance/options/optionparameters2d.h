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


#ifndef APPLICATIONS_FINANCE_OPTIONS_OPTIONPARAMETERS2D_H
#define APPLICATIONS_FINANCE_OPTIONS_OPTIONPARAMETERS2D_H 1


#include <applications/finance/options/optiontypesnd.h>

namespace lawa {

template < typename T, OptionTypenD OType>
struct OptionParameters2D
{

};

template <typename T>
struct OptionParameters2D<T,BasketPut> {

    OptionParameters2D(void);

    OptionParameters2D(T _strike, T _maturity, T weight1, T weight2, bool _earlyExercise);

    T strike;
    T maturity;
    T weight1, weight2;
    bool earlyExercise;
};

template <typename T>
struct OptionParameters2D<T,SumOfPuts> {

    OptionParameters2D(void);

    OptionParameters2D(T _strike1, T strike2, T _maturity, T weight1, T weight2, bool _earlyExercise);

    T strike1, strike2;
    T maturity;
    T weight1, weight2;
    bool earlyExercise;
};

}   // namespace lawa

#include <applications/finance/options/optionparameters2d.tcc>

#endif  // APPLICATIONS_FINANCE_OPTIONS_OPTIONPARAMETERS2D_H
