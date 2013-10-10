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

#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_ADAPTIVESEPARABLERHS_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_ADAPTIVESEPARABLERHS_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3=RHSIntegral_x2, typename RHSIntegral_x4=RHSIntegral_x2,
          typename RHSIntegral_x5=RHSIntegral_x2>
class AdaptiveSeparableRhs
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;

    public:

    AdaptiveSeparableRhs(const RHSIntegral_x1 &_integral_x1, Coefficients<Lexicographical,T,Index1D> &_f_data_x1,
                         const RHSIntegral_x2 &_integral_x2, Coefficients<Lexicographical,T,Index1D> &_f_data_x2);

    AdaptiveSeparableRhs(const RHSIntegral_x1 &_integral_x1, Coefficients<Lexicographical,T,Index1D> &_f_data_x1,
                         const RHSIntegral_x2 &_integral_x2, Coefficients<Lexicographical,T,Index1D> &_f_data_x2,
                         const RHSIntegral_x3 &_integral_x3, Coefficients<Lexicographical,T,Index1D> &_f_data_x3);


    T
    operator()(const Index2D &index2d);

    T
    operator()(const Index3D &index3d);

    void
    clear();

    const RHSIntegral_x1                    &integral_x1;
    Coefficients<Lexicographical,T,Index1D> &f_data_x1;
    const RHSIntegral_x2                    &integral_x2;
    Coefficients<Lexicographical,T,Index1D> &f_data_x2;
    const RHSIntegral_x3                    &integral_x3;
    Coefficients<Lexicographical,T,Index1D> &f_data_x3;
    const RHSIntegral_x4                    &integral_x4;
    Coefficients<Lexicographical,T,Index1D> &f_data_x4;
    const RHSIntegral_x5                    &integral_x5;
    Coefficients<Lexicographical,T,Index1D> &f_data_x5;
};

}   // namespace lawa

#include <lawa/methods/adaptive/righthandsides/adaptiveseparablerhs.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_ADAPTIVESEPARABLERHS_H
