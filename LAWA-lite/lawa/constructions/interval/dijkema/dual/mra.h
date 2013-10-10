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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/interval/refinementmatrix.h>

namespace lawa {

using namespace flens;

template <typename _T>
class MRA<_T,Dual,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Dual;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef BasisFunction<T,Dual,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Dual,Interval,Dijkema> BSplineType;

        MRA(int d, int d_, int j=-1);

        ~MRA();

        // cardinalities of whole, left, inner, right index sets.
        int
        cardI_(int j) const;

        int
        cardI_L(int j=0) const;

        int
        cardI_I(int j) const;

        int
        cardI_R(int j=0) const;

        // ranges of whole left, inner, right index sets.
        Range<int>
        rangeI_(int j) const;

        Range<int>
        rangeI_L(int j=0) const;

        Range<int>
        rangeI_I(int j) const;

        Range<int>
        rangeI_R(int j) const;

        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const int d, d_, mu;   // mu = mu(d) = d&1.

    private:
        const int l1, l2;        // support of phi  = [ l1, l2 ] (real line).
        const int l1_, l2_;      // support of phi  = [ l1, l2 ] (real line).

    public:
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

        BSpline<T,Dual,Interval,Dijkema> phi_;
        GeMatrix<FullStorage<T,cxxblas::ColMajor> > R_Left, R_Right;
        RefinementMatrix<T,Interval,Dijkema> M0_;
        BSpline<T,Dual,R,CDF> phi_R;

    private:
        void
        _calcM0_();

    public: // FIXME: "public: " TO BE ELIMINATED !!!!!!!!!!!!!!
        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/dual/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H

