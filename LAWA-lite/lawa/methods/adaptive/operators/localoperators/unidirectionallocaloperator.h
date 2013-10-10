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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALLOCALOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALLOCALOPERATOR_H 1

#include <cstring>
#include <lawa/flensforlawa.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/operators/localoperators/localoperator1d.h>

namespace lawa {

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
struct UniDirectionalLocalOperator
{
    typedef typename LocalOperator1D::T T;

    typedef NotCoordXIndex                                                 notCoordXIndex;

    typedef typename LocalOperator1D::TrialWaveletBasis                    TrialBasis_CoordX;
    typedef typename LocalOperator1D::TestWaveletBasis                     TestBasis_CoordX;

    typedef typename TreeCoefficients1D<T>::const_by_level_it              const_by_level_it;
    typedef typename TreeCoefficients1D<T>::by_level_it                    by_level_it;

    typedef AlignedCoefficients<T,Index,NotCoordXIndex,Index1D,NotCoordX>  NotCoordXAlignedCoefficients;

    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    typedef typename IndexSet<Index1D>::const_iterator                       const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator   const_coeff_it;

    UniDirectionalLocalOperator(LocalOperator1D &_localOperator1D, T _factor=1.);

    void
    setParameters(int _J, size_t _hashTableLargeLength, size_t _hashTableSmallLength);

    void
    eval(const Coefficients<Lexicographical,T,Index> &v,
         Coefficients<Lexicographical,T,Index> &IAIv);

    void
    eval(const Coefficients<Lexicographical,T,Index> &v,
         Coefficients<Lexicographical,T,Index> &IAIv, const char* evalType);

    void
    nonTreeEval(const Index1D &coordX_col_index, const NotCoordXIndex &notcoordX_col_index,
                T col_val, IndexSet<Index1D> &row_indices1d,
                Coefficients<Lexicographical,T,Index> &Av, T eps=0.);

    LocalOperator1D          &localOperator1D;
    const TrialBasis_CoordX  &trialBasis_CoordX;
    const TestBasis_CoordX   &testBasis_CoordX;
    int                      J;
    size_t                   hashTableLargeLength;
    size_t                   hashTableSmallLength;

    T                        factor;

    Split<Index,Index1D,NotCoordXIndex,CoordX> split;
    Join<Index,Index1D,NotCoordXIndex,CoordX>  join;

};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/localoperators/unidirectionallocaloperator.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALLOCALOPERATOR_H
