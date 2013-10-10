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


#ifndef  LAWA_RIGHTHANDSIDES_RHSWITHPEAKS1D_H
#define  LAWA_RIGHTHANDSIDES_RHSWITHPEAKS1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/functiontypes/function.h>

namespace lawa {

/*
 * Computes <f, psi_{j,k}> + sum_{i=1}^n deltas(i,2) * psi_{j,k}(deltas(i,1))
 * _with_singular_part:  indicates whether sum_{i=1}^n deltas(i,2) * psi_{j,k}(deltas(i,1))
 *                       is computed or not
 * _with_smooth_part:    indicates whether  <f, psi_{j,k}> is computed or not
 */

template <typename T, typename Basis>
class RHSWithPeaks1D
{
    public:
        RHSWithPeaks1D(const Basis &basis, Function<T> _f,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas,
                       int order, bool _with_singular_part=true, bool _with_smooth_part=true);

        T
        operator()(XType xtype, int j, long k) const;

        T
        operator()(const Index1D &lambda) const;

    private:
        const Basis &basis;
        Function<T> f;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    &deltas;
        bool with_singular_part;
        bool with_smooth_part;
        IntegralF<Gauss, Basis> integralf;

};

template <typename T>
class RHSWithPeaks1D_WO_XBSpline
{
    public:

        RHSWithPeaks1D_WO_XBSpline(const Wavelet<T,Primal,R,CDF> &_psi,
                                   T (*_f)(T), const DenseVector<Array<T> > &_f_singularPoints,
                                   flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &_deltas,
                                   T _left_bound, T _right_bound, T _h, int _order,
                                   bool _with_singular_part=true, bool _with_smooth_part=true);
        T
        truncated_f(T x) const;

        T
        operator()(const Index1D &lambda) const;

        T
        operator()(int j, long k, T a, T b) const;

    private:

        const Wavelet<T,Primal,R,CDF> &psi;
        T (*f)(T);
        const DenseVector<Array<T> > f_singularPoints;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &deltas;
        T left_bound, right_bound;
        T h;    //interval length for quadrature
        int order;
        bool with_singular_part;
        bool with_smooth_part;

        DenseVector<Array<T> >       f_singularPoints_interval;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;

};


} // namespace lawa

#include <lawa/righthandsides/rhswithpeaks1d.tcc>

#endif //  LAWA_RIGHTHANDSIDES_RHSWITHPEAKS1D_H

