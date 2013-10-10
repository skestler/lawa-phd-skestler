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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H 1

#include <lawa/settings/enum.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/refinementmatrix.h>
#include <lawa/constructions/periodic/primal/mra.h>
#include <lawa/constructions/periodic/dual/mra.h>
#include <lawa/constructions/periodic/primal/wavelet.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Primal,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Primal,Periodic,CDF> BasisFunctionType;
        typedef BSpline<T,Primal,Periodic,CDF> BSplineType;
        typedef Wavelet<T,Primal,Periodic,CDF> WaveletType;

        Basis(int _d, int _d_, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        const BasisFunctionType &
        generator(XType xtype) const;
        
        int
        cardJ(int j) const;

        Range<int>
        rangeJ(int j) const;

        const int d, d_, j0;
        MRA<T,Primal,Periodic,CDF> mra;
        MRA<T,Dual,Periodic,CDF> mra_;
        Wavelet<T,Primal,Periodic,CDF> psi;
        RefinementMatrix<T,Periodic,CDF> M1;
        
    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H

