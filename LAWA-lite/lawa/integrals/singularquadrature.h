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

#ifndef LAWA_INTEGRALS_SINGULARQUADRATURE_H
#define LAWA_INTEGRALS_SINGULARQUADRATURE_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {


template <typename SingularIntegral>
class SingularQuadrature
{
    typedef flens::DenseVector<flens::Array<long double> >                        DenseVector;

    public:
        SingularQuadrature(const SingularIntegral &_singularintegral);

        void
        setLegendreOrder(int order_eta);

        void
        setParameters(int order, int n, double sigma, double mu, double omega);

        const long double
        operator()(long double a1, long double b1, long double a2, long double b2,
                   long double eps);

        const SingularIntegral &singularintegral;

    private:
        long double
        _integrate_singular_diagonal(long double a1, long double b1, long double a2, long double b2,
                                     long double eps);

        long double
        _integrate_singular_corner(long double a1, long double b1, long double a2, long double b2,
                                   long double eps);

        long double
        _integrate_nonsingular(long double a1, long double b1, long double a2, long double b2,
                               long double eps);

        static void
        _legendre(int order);

        static void
        _hp_composite_legendre(int n, double sigma, double mu);

        int _order_eta;      //Legendre order for tensor product rule
        int _order;
        int _n;              //parameters for composite variable order Gauss-Legendre quadrature
        double _sigma;
        double _mu;
        double _omega;
        static int _precalculated_order;
        static int _precalculated_n;
        static double   _precalculated_sigma;
        static double   _precalculated_mu;
        static flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > _legendreknots;
        static flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > _legendreweights;
        static flens::DenseVector<flens::Array<int> >                              _hp_legendrenumofpoints;
        static flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > _hp_legendreknots;
        static flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > _hp_legendreweights;
};



} // namespace lawa

#include <lawa/integrals/singularquadrature.tcc>

#endif // LAWA_INTEGRALS_SINGULARQUADRATURE_H

