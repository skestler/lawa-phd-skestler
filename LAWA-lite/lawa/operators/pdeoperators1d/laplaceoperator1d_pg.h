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


#ifndef LAWA_OPERATORS_PDEOPERATORS1D_LAPLACEOPERATOR1D_PG_H
#define LAWA_OPERATORS_PDEOPERATORS1D_LAPLACEOPERATOR1D_PG_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>

namespace lawa {

/* Laplace OPERATOR
 *
 *    a(v,u) =  Integral(v_x * u_x)
 *
 */
template <typename T, typename TrialBasis, typename TestBasis>
class LaplaceOperator1D_PG{
    
    public:
        
        const TrialBasis& trialbasis;
        const TestBasis& testbasis;
      
        LaplaceOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis);

        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
        
    private:

        Integral<Gauss, TestBasis, TrialBasis> integral;

};

} // namespace lawa

#include <lawa/operators/pdeoperators1d/laplaceoperator1d_pg.tcc>

#endif  //LAWA_OPERATORS_PDEOPERATORS1D_LAPLACEOPERATOR1D_PG_H 1

