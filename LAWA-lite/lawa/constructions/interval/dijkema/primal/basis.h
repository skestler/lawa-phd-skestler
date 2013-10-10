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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Primal,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef Basis<T,Primal,Interval,Dijkema>         RefinementBasis;
        typedef BasisFunction<T,Primal,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Dijkema> BSplineType;
        typedef Wavelet<T,Primal,Interval,Dijkema> WaveletType;

        Basis(int _d, int _d_, int j=-1);
        
        ~Basis();

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
        cardJ(int j) const;

        int
        cardJL(int j=-1) const;

        int
        cardJI(int j) const;

        int
        cardJR(int j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const Range<int>
        rangeJ(int j) const;

        const Range<int>
        rangeJL(int j=-1) const;

        const Range<int>
        rangeJI(int j) const;

        const Range<int>
        rangeJR(int j=-1) const;

        /// Returns the range of indicicated functions wavelets from SecondBasis or
        /// SecondRefinemementBasiswhose supports intersect the support
        /// of a refinement B-spline with level j_bspline and translation index k_bspline
        /// from the current RefinementBasis. This is required for tree-based algorithms.
        /// The returned level j_wavelet or j_bspline2 is chosen s.t. there is "no scale difference"
        /// i.e., if we refine both functions the corresponding refinements should
        /// live on the same scale.
        template <typename SecondRefinementBasis>
            void
            getWaveletNeighborsForBSpline(int j_bspline, long k_bspline,
                                          const SecondRefinementBasis &secondbasis,
                                          int &j_wavelet, long &k_wavelet_first,
                                          long &k_wavelet_last) const;

        template <typename SecondRefinementBasis>
            void
            getBSplineNeighborsForBSpline(int j_bspline1, long k_bspline1,
                                          const SecondRefinementBasis &secondrefinementbasis,
                                          int &j_bspline2,
                                          long &k_bspline2_first, long &k_bspline2_last) const;


        /// Returns the range of indicated functions from SecondBasis whose supports
        /// intersect the support of a given (multi-)scaling with level j_scaling and translation index
        /// k_scaling from the current Basis. This is required for tree-based algorithms.
        /// The returned level is chosen s.t. the corresponding refinements live
        /// on the same scale.
        template <typename SecondBasis>
            void
            getScalingNeighborsForScaling(int j_scaling1, long k_scaling1,
                                          const SecondBasis &secondbasis,
                                          int &j_scaling2, long &k_scaling_first,
                                          long &k_scaling_last) const;

        template <typename SecondBasis>
            void
            getWaveletNeighborsForScaling(int j_scaling1, long k_scaling1,
                                          const SecondBasis &secondbasis,
                                          int &j_wavelet, long &k_wavelet_first,
                                          long &k_wavelet_last) const;

        /// Returns the range of indicated functions from SecondBasis and SecondRefinementBasis
        /// whose supports intersect the support of a given wavelet with level j_wavelet and
        /// translation index k_wavelet from the current Basis. This is required for tree-based algorithms.
        /// The returned level of the functions is chosen s.t. there is "no scale difference", i.e.,
        /// the corresponding refinements should live on the same scale.
        template <typename SecondBasis>
            void
            getBSplineNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                          const SecondBasis &secondrefinementbasis,
                                          int &j_bspline, long &k_bspline_first,
                                          long &k_bspline_last) const;

        template <typename SecondBasis>
            void
            getScalingNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                          const SecondBasis &secondbasis,
                                          int &j_scaling, long &k_scaling_first,
                                          long &k_scaling_last) const;

        template <typename SecondBasis>
            void
            getWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                          const SecondBasis &secondbasis,
                                          int &j_wavelet2, long &k_wavelet_first,
                                          long &k_wavelet_last) const;

        template <typename SecondBasis>
            void
            getLowerWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                               const SecondBasis &secondbasis,
                                               int &j_wavelet2, long &k_wavelet_first,
                                               long &k_wavelet_last) const;

        template <typename SecondBasis>
            void
            getHigherWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                               const SecondBasis &secondbasis,
                                               int &j_wavelet2, long &k_wavelet_first,
                                               long &k_wavelet_last) const;

        class LaplaceOperator1D {
            public:
                LaplaceOperator1D(int _d,
                                  const Basis<_T,Primal,Interval,Dijkema> &_refinementbasis);

                T
                operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2);

            private:
                int d;
                const Basis<_T,Primal,Interval,Dijkema> &refinementbasis;

                DenseVector<Array<long double> > outer_values;
                DenseVector<Array<long double> > inner_values;
        };

        class IdentityOperator1D {
            public:
                IdentityOperator1D(int _d,
                                   const Basis<_T,Primal,Interval,Dijkema> &_refinementbasis);

                T
                operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2);

            private:
                int d;
                const Basis<_T,Primal,Interval,Dijkema> &refinementbasis;

                DenseVector<Array<long double> > outer_values;
                DenseVector<Array<long double> > inner_values;
        };

        LaplaceOperator1D LaplaceOp1D;
        IdentityOperator1D IdentityOp1D;

        MRA<T,Primal,Interval,Dijkema> mra;
        MRA<T,Dual,Interval,Dijkema>  mra_;

        RefinementMatrix<T,Interval,Dijkema> M1;

        const int d, d_, mu;   // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

    private:
        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

        friend class Wavelet<T,Primal,Interval,Dijkema>;

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
        Wavelet<T,Primal,Interval,Dijkema> psi;

        Basis<T,Primal,Interval,Dijkema> &refinementbasis;
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H

