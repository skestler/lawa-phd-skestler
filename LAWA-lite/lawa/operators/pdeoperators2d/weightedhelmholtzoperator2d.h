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

#ifndef LAWA_OPERATORS_PDEOPERATORS2D_WEIGHTEDHELMHOLTZOPERATOR2D_H
#define LAWA_OPERATORS_PDEOPERATORS2D_WEIGHTEDHELMHOLTZOPERATOR2D_H 1

#include <lawa/functiontypes/function.h>
#include <lawa/integrals/integral.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/uniformoperator2d.h>

namespace lawa {

/* Weighted Laplace OPERATOR 2D 
 *
 *    a(v,u) =  	Integral(w1 * v1_x * u1_x) * Integral(w2 * v2 * u2) 
 *			  + 	Integral(w1 * v1 * u1) 	   * Integral(w2 * v2_y * u2_y) 
 *			  + c * Integral(w1 * v1 * u1) 	   * Integral(w2 * v2 * u2);
 *
 */
template <typename T, typename Basis>
class WeightedHelmholtzOperator2D : public UniformOperator2D<T> {
	
	public:
  		const Basis& basis;
        const T c;
        
        WeightedHelmholtzOperator2D(const Basis& _basis, const T _c, Function<T> weightFct_x, Function<T> weightFct_y);
        
        T
    	operator()(XType row_xtype_x, int j1_x, long k1_x,
                   XType row_xtype_y, int j1_y, long k1_y,
                   XType col_xtype_x, int j2_x, long k2_x,
                   XType col_xtype_y, int j2_y, long k2_y) const;
    	
        T
        operator()(const Index2D &row_index, const Index2D &col_index) const;
        
        T
        operator()(const Index2D &row_index, const Index2D &col_index);
    
    private:

        typedef typename Basis::FirstBasisType Basis_x;
        typedef typename Basis::SecondBasisType Basis_y;
        
        Function<T> W_x;
        Function<T> W_y;
        
        IntegralF<Gauss, Basis_x, Basis_x>   integral_x;
        IntegralF<Gauss, Basis_y, Basis_y>   integral_y;    
};

} // namespace lawa

#include <lawa/operators/pdeoperators2d/weightedhelmholtzoperator2d.tcc>

#endif // LAWA_OPERATORS_PDEOPERATORS2D_WEIGHTEDHELMHOLTZOPERATOR2D_H
