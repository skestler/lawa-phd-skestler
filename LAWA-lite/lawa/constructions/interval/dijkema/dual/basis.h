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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Dual,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Dual;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef BasisFunction<T,Dual,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Dual,Interval,Dijkema> BSplineType;
        typedef Wavelet<T,Dual,Interval,Dijkema> WaveletType;

        Basis(int _d, int _d_, int j=-1);

        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const BasisFunctionType &
        generator(XType xtype) const;

        // cardinalities of whole, left, inner, right index sets (primal).
        int
        cardJ_(int j) const;

        int
        cardJ_L(int j=-1) const;

        int
        cardJ_I(int j) const;

        int
        cardJ_R(int j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const Range<int>
        rangeJ_(int j) const;

        const Range<int>
        rangeJ_L(int j=-1) const;

        const Range<int>
        rangeJ_I(int j) const;

        const Range<int>
        rangeJ_R(int j=-1) const;

        MRA<T,Primal,Interval,Dijkema> mra;
        MRA<T,Dual,Interval,Dijkema>  mra_;

        RefinementMatrix<T,Interval,Dijkema> M1_;

        const int d, d_, mu;   // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level. 

    private:
        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

    public:
        Wavelet<T,Dual,Interval,Dijkema> psi_;
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/dual/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H

