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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/primal/bspline.h>
#include <lawa/constructions/periodic/refinementmatrix.h>

namespace lawa {

template <typename _T>
class MRA<_T,Primal,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Primal,Periodic,CDF> BasisFunctionType;
        typedef BSpline<T,Primal,Periodic,CDF> BSplineType;

        MRA(int _d, int _d_, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        int
        cardI(int j) const;

        Range<int>
        rangeI(int j) const;

        const int d, d_, j0;
        BSpline<T,Primal,Periodic,CDF> phi;
        RefinementMatrix<T,Periodic,CDF> M0;
        
    protected:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H

