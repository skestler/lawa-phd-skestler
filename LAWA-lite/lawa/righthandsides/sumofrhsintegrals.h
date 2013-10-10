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

#ifndef LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS_H
#define LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS_H 1

namespace lawa {

template <typename T, typename Index, typename RHSIntegral1, typename RHSIntegral2>
class SumOfTwoRHSIntegrals
{
private:
    const RHSIntegral1 &rhs1;
    const RHSIntegral2 &rhs2;

public:
    SumOfTwoRHSIntegrals(const RHSIntegral1 &rhs1, const RHSIntegral2 &rhs2);

    T
    operator()(const Index &index) const;
    
    T
    operator()(XType xtype_x, int j_x, long k_x,
               XType xtype_y, int j_y, long k_y) const;
};


template <typename T, typename Index, typename RHSIntegral>
class SumOfThreeRHSIntegrals
{
private:
    const RHSIntegral &rhs1;
    const RHSIntegral &rhs2;
    const RHSIntegral &rhs3;

public:
    SumOfThreeRHSIntegrals(const RHSIntegral &rhs1, const RHSIntegral &rhs2, const RHSIntegral &rhs3);

    T
    operator()(const Index &index) const;
};


}    //namespace lawa

#include <lawa/righthandsides/sumofrhsintegrals.tcc>

#endif  //LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS2D_H
