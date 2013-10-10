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

#ifndef LAWA_INTEGRALS_QUADRATURE_H
#define LAWA_INTEGRALS_QUADRATURE_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <QuadratureType Quad, typename Integral>
class Quadrature
{
};

template <typename Integral>
class Quadrature<Gauss,Integral>
{
    public:
        typedef typename Integral::T T;

        Quadrature(const Integral &_integral);

        const T
        operator()(T a, T b) const;

        int
        order() const;

        void
        setOrder(int order);

        const Integral &integral;
    private:
        static void
        _legendre(int order);

        int _order;
        static int _precalculatedOrder;
        static flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        static flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

//-----------------------------------------------------------------------------

template <typename Integral>
class Quadrature<Trapezoidal,Integral>
{
    public:
        typedef typename Integral::T T;
        
        Quadrature(const Integral &_integral);

        const T
        operator()(T a, T b) const;

        int
        n() const;

        void
        setN(int n);

        const Integral &integral;
    private:
        int _n;
};

//-----------------------------------------------------------------------------

template <typename Integral>
class Quadrature<ExpWeighted,Integral>
{
    public:
        typedef typename Integral::T T;

        Quadrature(const Integral &_integral);

        const T
        operator()(T a, T b) const;

        const Integral &integral;

    private:
        const T _eta;
        T _max_polynomialorder;
};

} // namespace lawa

#include <lawa/integrals/quadrature.tcc>

#endif // LAWA_QUADRATURE_H

