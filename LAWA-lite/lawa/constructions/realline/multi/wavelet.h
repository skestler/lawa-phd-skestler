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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>

namespace lawa {

template <typename _T>
class Wavelet<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
    
        Wavelet(int _d);
        
        Wavelet(const Basis<T,Orthogonal,R,Multi> &basis);
        
        virtual
        ~Wavelet();
        
        T
        operator()(T x, int j, long k, unsigned short deriv) const;
        
        Support<T>
        support(int j, long k) const;
        
        Support<T>
        max_support() const;

        DenseVector<Array<T> >
        singularSupport(int j, long k) const;
        
        T
        tic(int j) const;

        DenseVector<Array<long double> > *
        getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const;

        int
        getRefinementLevel(int j) const;

        const int d;
        const int vanishingMoments;
        unsigned int _numSplines;

//    private:      // should be private one fine day

        typedef T (*Evaluator)(T x, unsigned short deriv);

        void
        _initialize(int d);
        
        long
        _shift(long k) const;
        
        int
        _type(long k) const;
    
        Evaluator *_evaluator;
        Support<T> *_support;
        DenseVector<Array<T> > *_singularSupport;

        Support<T> _max_support;

        DenseVector<Array<long double> > *_refCoeffs;
        long double                      *_H1SemiNorms;
        long                             *_offsets;
        T                                _initialticsize;
        int                              _addRefinementLevel;    //B-splines for refinement are needed on higher levels
        int                              _shiftFactor;           //Needed since we have multiple B-spline generators for refinement.
    
};
    
} // namespace lawa

#include <lawa/constructions/realline/multi/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H
