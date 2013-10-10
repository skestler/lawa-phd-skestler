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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDBASKETPUTOPTION2D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDBASKETPUTOPTION2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {


template<typename T>
struct TruncatedBasketPutOption2D
{
    static T   left_x1;
    static T   right_x1;
    static T   left_x2;
    static T   right_x2;
    static DenseVector<Array<T> > singPts_x1;
    static DenseVector<Array<T> > singPts_x2;

    static Option2D<T,BasketPut>  basketputoption;

    static T   u11, u21, u12, u22;

    static int type;
    static T   truncWidth;
    static T   damping_c;

    static DenseVector<Array<T> > critical_line_x1;
    static bool                   critical_above_x1;
    static DenseVector<Array<T> > critical_line_x2;
    static bool                   critical_above_x2;


    static void
    setOption(const Option2D<T,BasketPut> &_basketputoption);

    static void
    setTransformation(T _u11, T _u21, T _u12, T _u22);

    static void
    setTruncation(T _left_x1, T _right_x1, T _left_x2, T _right_x2, int _type, T _truncWidth, T _damping_c);

    static void
    setCriticalLine_x1(T _critical_line_x1, bool critical_above_x1);

    static void
    setCriticalLine_x2(T _critical_line_x2, bool critical_above_x2);

    static bool
    isCritical(T a1, T b1, T a2, T b2);

    static T
    payoff(T x, T y);
};

}   // namespace lawa

#include <applications/finance/initialconditions/truncatedbasketputoption2d.tcc>

#endif // APPLICATIONS_FINANCE_INITIALCONDITIONS_TRUNCATEDBASKETPUTOPTION2D_H
