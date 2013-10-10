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


#ifndef LAWA_OPERATORS_SPACETIMEOPERATORS_SPACETIMEINITIALCONDITION1D_H
#define LAWA_OPERATORS_SPACETIMEOPERATORS_SPACETIMEINITIALCONDITION1D_H 1


#include <lawa/integrals/integral.h>

namespace lawa {

/* Space-Time Initial Condition
 *
 *      a(v,u) = u1(0) * Integral(v2 * u2)
 * 
 */
template <typename T, typename Basis>
class SpaceTimeInitialCondition1D{

    private:

        typedef typename Basis::FirstBasisType Basis_t;
        typedef typename Basis::SecondBasisType Basis_x;

        Integral<Gauss, Basis_x, Basis_x>   integral_x;

    public:
        
        const Basis& basis;
                
        SpaceTimeInitialCondition1D(const Basis& _basis);

        T                                                           // returns a(v,u)
        operator()(XType row_xtype_x, int j1_x, long k1_x,
                   XType col_xtype_t, int j2_t, long k2_t,
                   XType col_xtype_x, int j2_x, long k2_x) const;

        T
        operator()(const Index1D &row_index, const Index2D &col_index) const;
};


} // namespace lawa

#include <lawa/operators/spacetimeoperators/spacetimeinitialcondition1d.tcc>

#endif // LAWA_OPERATORS_SPACETIMEOPERATORS_SPACETIMEINITIALCONDITION1D_H

