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


#ifndef LAWA_OPERATORS_PDEOPERATORS1D_HELMHOLTZOPERATOR1D_H
#define LAWA_OPERATORS_PDEOPERATORS1D_HELMHOLTZOPERATOR1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>
#include <lawa/settings/enum.h>

namespace lawa {
    
/* HELMHOLTZ OPERATOR 
 *
 *    a(v,u) =  Integral(v_x * u_x) + c * Integral(v * u)
 *
 */
template <typename T, typename Basis>
class HelmholtzOperator1D
{
    public:

        const Basis& basis;
        const T c;

        HelmholtzOperator1D(const Basis& _basis, const T _c);
        HelmholtzOperator1D(const HelmholtzOperator1D<T,Basis> &a);
    
        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;

        T
        operator()(T time,
                   XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(T time, const Index1D &row_index, const Index1D &col_index) const;

    private:

        Integral<Gauss, Basis, Basis> integral;
};

} // namespace lawa

#include <lawa/operators/pdeoperators1d/helmholtzoperator1d.tcc>

#endif // LAWA_OPERATORS_PDEOPERATORS1D_HELMHOLTZOPERATOR1D_H

