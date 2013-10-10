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

#ifndef LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_WAVELET_H
#define LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/wavelet.h>
#include <cassert>


namespace lawa {

template <typename _T>
class Wavelet<_T,Primal,RPlus,SparseMulti>
    : public BasisFunction<_T,Primal,RPlus,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = RPlus;
        static const Construction Cons = SparseMulti;
            
        Wavelet(const Basis<T,Primal,RPlus,SparseMulti> &basis);
        
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
        
        const Basis<T,Primal,RPlus,SparseMulti> &basis;
        const int d;
        const int vanishingMoments;
        int _numSplines;
        Support<T> _max_support;
};
    
} // namespace lawa

#include <lawa/constructions/rplus/sparsemulti/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_WAVELET_H
