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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Wavelet<_T,Orthogonal,Interval,Multi>
    : public BasisFunction<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;
            
        Wavelet(const Basis<T,Orthogonal,Interval,Multi> &basis);
        
        virtual
        ~Wavelet();
        
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

        T
        getL2Norm(int j, long k) const;

        T
        getH1SemiNorm(int j, long k) const;

        const Basis<T,Orthogonal,Interval,Multi> &basis;
        const int d;
        const int vanishingMoments;
        T initialticsize;
};
    
} // namespace lawa

#include <lawa/constructions/interval/multi/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_H
