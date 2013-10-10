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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_MRA_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/interval/refinementmatrix.h>

namespace lawa {

using namespace flens;

template <typename _T>
class MRA<_T,Primal,Interval,Primbs>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Primbs;

        typedef BasisFunction<T,Primal,Interval,Primbs> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Primbs> BSplineType;

        MRA(int d, int j=-1);

        ~MRA();

        T
        operator()(const T x, const int j, const int k) const;

        // cardinalities of whole, left, inner, right index sets.
        int
        cardI(int j) const;

        int
        cardIL(int j=0) const;

        int
        cardII(int j) const;

        int
        cardIR(int j=0) const;

        // ranges of whole left, inner, right index sets.
        Range<int>
        rangeI(int j) const;

        Range<int>
        rangeIL(int j=0) const;

        Range<int>
        rangeII(int j) const;

        Range<int>
        rangeIR(int j) const;

        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const int d, mu;       // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

        BSpline<T,Primal,R,CDF> phiR;
        RefinementMatrix<T,Interval,Primbs> M0;

        const int l1, l2;      // support of phi  = [ l1, l2 ] (real line).
        
    private:
        void
        _calcM0();

        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

        friend class BSpline<T,Primal,Interval,Dijkema>;

        DenseVector<Array<long double> > *_leftRefCoeffs,
                                         *_innerRefCoeffs,
                                         *_rightRefCoeffs;

        long double *_leftL2Norms,  *_leftH1SemiNorms,
                    *_innerL2Norms, *_innerH1SemiNorms,
                    *_rightL2Norms, *_rightH1SemiNorms;

        long *_leftOffsets,
             *_innerOffsets,
             *_rightOffsets;

    public:
        BSpline<T,Primal,Interval,Primbs> phi;

};

} // namespace lawa

#include <lawa/constructions/interval/primbs/primal/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_MRA_H

