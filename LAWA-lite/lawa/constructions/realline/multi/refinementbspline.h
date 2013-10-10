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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,R,MultiRefinement>
    : public BasisFunction<_T,Orthogonal,R,MultiRefinement>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = MultiRefinement;

        BSpline(const int _d);

        virtual
        ~BSpline();

        T
        operator()(T x, int j, long k, unsigned short deriv) const;

        Support<T>
        support(int j, long k) const;

        DenseVector<Array<T> >
        singularSupport(int j, long k) const;

        T
        tic(int j) const;

        DenseVector<Array<long double> > *
        getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const;

        int
        getRefinementLevel(int j) const;

        const unsigned int d;
        unsigned int _numSplines;

    //private:  // should be private one fine day...

        typedef T (*Evaluator)(T x, unsigned short deriv);

        long
        _shift(long k) const;

        int
        _type(long k) const;

        Evaluator                        *_evaluator;
        Support<T>                       *_support;
        DenseVector<Array<T> >           *_singularSupport;
        DenseVector<Array<long double> > *_refCoeffs;
        long                             *_offsets;
        T                                _initialticsize;
        int                              _shiftFactor;
};

} // namespace lawa

#include <lawa/constructions/realline/multi/refinementbspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBSPLINE_H
