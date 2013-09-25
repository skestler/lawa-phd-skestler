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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDPUTOPTION1D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDPUTOPTION1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {


template<typename T, OptionType1D OType>
struct TruncatedPutOption1D
{
    static T   left;
    static T   right;
    static int type;
    static T   truncWidth;

    static DenseVector<Array<T> > singPts;
    static Option1D<T,OType>  option;

    static void
    setOption(const Option1D<T,Put> &_option);

    static void
    setTruncation(T left, T right, int _type, T _truncWidth);

    static T
    g_trunc(T x);
};

}   // namespace lawa

#include <applications/finance/initialconditions/truncatedputoption1d.tcc>

#endif // APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDPUTOPTION1D_H
