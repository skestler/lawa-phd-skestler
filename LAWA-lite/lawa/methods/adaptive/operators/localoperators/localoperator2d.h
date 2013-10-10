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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_LOCALOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_LOCALOPERATOR2D_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/operators/localoperators/localoperator1d.h>

namespace lawa {

template <typename LocalOperator1, typename LocalOperator2>
struct LocalOperator2D {

    typedef typename LocalOperator1::T T;

    typedef typename LocalOperator1::TrialWaveletBasis                         TrialBasis_x1;
    typedef typename LocalOperator1::TestWaveletBasis                          TestBasis_x1;
    typedef typename LocalOperator2::TrialWaveletBasis                         TrialBasis_x2;
    typedef typename LocalOperator2::TestWaveletBasis                          TestBasis_x2;


    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef IndexSet<Index2D>::const_iterator                                  const_set2d_it;

    typedef typename Coefficients<Lexicographical,T,Index1D>::iterator         coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator   const_coeff2d_it;

    typedef typename TreeCoefficients1D<T>::const_by_level_it                  const_by_level_it;
    typedef typename TreeCoefficients1D<T>::by_level_it                        by_level_it;

    typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XOne>             XOneAlignedCoefficients;
    typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XTwo>             XTwoAlignedCoefficients;

    LocalOperator2D(LocalOperator1 &_localoperator1, LocalOperator2 &_localoperator2);

    void
    setJ(int J);

    void
    eval(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &AAv);

    void
    eval(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &AAv,
         T &time_intermediate1, T &time_intermediate2,
         T &time_IAv1, T &time_IAv2, T &time_LIv, T &time_UIv);

    void
    debug_eval(const Coefficients<Lexicographical,T,Index2D> &v,
                 Coefficients<Lexicographical,T,Index2D> &IAUIv,
                 Coefficients<Lexicographical,T,Index2D> &LIIAv,
                 const Coefficients<Lexicographical,T,Index2D> &IAv_ref,
                 const Coefficients<Lexicographical,T,Index2D> &LIIAv_ref,
                 const Coefficients<Lexicographical,T,Index2D> &UIv_ref,
                 const Coefficients<Lexicographical,T,Index2D> &IAUIv_ref,
                 const Coefficients<Lexicographical,T,Index2D> &AAv_ref) /*const*/;

    void
    initializeIntermediateVectorIAv(const Coefficients<Lexicographical,T,Index2D> &v,
                                    const Coefficients<Lexicographical,T,Index2D> &LIIAv,
                                    Coefficients<Lexicographical,T,Index2D> &IAv) const;

    void
    initializeIntermediateVectorUIv(const Coefficients<Lexicographical,T,Index2D> &v,
                                    const Coefficients<Lexicographical,T,Index2D> &IAUIv,
                                    Coefficients<Lexicographical,T,Index1D> &Pe1_UIv) const;

    void
    evalIA(const Coefficients<Lexicographical,T,Index2D> &z,
           Coefficients<Lexicographical,T,Index2D> &IAz) /*const*/;

    void
    evalLI(const Coefficients<Lexicographical,T,Index2D> &z,
              Coefficients<Lexicographical,T,Index2D> &LIz) /*const*/;

    void
    evalUI(const Coefficients<Lexicographical,T,Index2D> &z,
           const Coefficients<Lexicographical,T,Index1D> &Pe1_UIz,
           Coefficients<Lexicographical,T,Index2D> &UIz) /*const*/;


    LocalOperator1          &localoperator1;
    LocalOperator2          &localoperator2;
    const TrialBasis_x1     &trialBasis_x1;
    const TestBasis_x1      &testBasis_x1;
    const TrialBasis_x2     &trialBasis_x2;
    const TestBasis_x2      &testBasis_x2;
    int                     J;
    size_t                  hashTableLargeLength;
    size_t                  hashTableSmallLength;

};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/localoperators/localoperator2d.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_LOCALOPERATOR2D_H
