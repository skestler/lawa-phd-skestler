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

#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H 1

#include <lawa/functiontypes/separablefunction2d.h>
#include <lawa/integrals/integrals.h>
#include <lawa/righthandsides/rhs2d.h>

/*
 * Computes ( <fx, psix_{jx,kx}> + sum_{i=1}^n deltasx(i,2) * psix_{j,k}(deltasx(i,1)) )
 *         +( <fy, psiy_{jy,ky}> + sum_{i=1}^m deltasy(i,2) * psiy_{j,k}(deltasy(i,1)) )
 */



namespace lawa {

template<typename T, typename Basis2D>
class SeparableRHS2D : public Rhs2D<T>
{
    public:
        SeparableRHS2D(const Basis2D& _basis, const SeparableFunction2D<T>& _F,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y,
                       int order);

        T
        operator()(XType xtype_x, int j_x, long k_x,
                   XType xtype_y, int j_y, long k_y) const;

        T
        operator()(const Index2D &index) const;


    private:
        const Basis2D& basis;
        const SeparableFunction2D<T>& F;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y;
        
        typedef typename Basis2D::FirstBasisType Basis_x;
        typedef typename Basis2D::SecondBasisType Basis_y;

        IntegralF<Gauss, Basis_x, Basis_x>  integralf_x;
        IntegralF<Gauss, Basis_y, Basis_y>  integralf_y;

};

} // namespace lawa

#include <lawa/righthandsides/separablerhs2d.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H

