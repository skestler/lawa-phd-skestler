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

#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_THETATIMESTEPSEPARABLERHS_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_THETATIMESTEPSEPARABLERHS_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
class ThetaTimeStepSeparableRHS {

    typedef typename Coefficients<Lexicographical,T,Index>::iterator           coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator     const_coeff_it;

    public:

    ThetaTimeStepSeparableRHS(Function<T> &_fct_t, SpatialRHS &_F_x,
                              ThetaTimeStepLocalOperator &_thetaTimeStepLocalOperator);

    void
    setThetaTimeStepParameters(T _theta,T _timestep, T _discrete_timepoint,
                               const Coefficients<Lexicographical,T,Index> &_u_k);

    void
    initializePropagation(const Coefficients<Lexicographical,T,Index> &f);

    T
    operator()(T t, const Index &index);

    T
    operator()(const Index &index);

    Function<T>                             &fct_t;
    SpatialRHS                              &F_x;
    ThetaTimeStepLocalOperator              &thetaTimeStepLocalOperator;
    Coefficients<Lexicographical,T,Index>   u_k;
    Coefficients<Lexicographical,T,Index>   propagated_u_k;
    T                                       discrete_timepoint;
    T                                       theta, timestep;

};

}   // namespace lawa

#include <lawa/methods/adaptive/righthandsides/thetatimestepseparablerhs.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_THETATIMESTEPSEPARABLERHS_H
