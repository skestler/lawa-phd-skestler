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


#ifndef APPLICATIONS_FINANCE_OPTIONS_OPTION_H
#define APPLICATIONS_FINANCE_OPTIONS_OPTION_H 1

#include <applications/finance/options/optiontypes1d.h>
#include <applications/finance/options/optiontypesnd.h>

namespace lawa {

template <typename T,OptionType1D Type>
struct Option1D
{

};

template <typename T,OptionTypenD Type>
struct Option2D
{

};

template <typename T,OptionTypenD Type>
struct Option3D
{

};

}   // namespace lawa

#endif  // APPLICATIONS_FINANCE_OPTIONS_OPTION_H
