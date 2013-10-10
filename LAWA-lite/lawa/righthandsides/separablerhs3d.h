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

#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H 1

#include <lawa/functiontypes/separablefunction3d.h>
#include <lawa/integrals/integrals.h>

namespace lawa {

template<typename T, typename Basis3D>
class SeparableRHS3D
{
    private:
        const Basis3D& basis;
        const SeparableFunction3D<T>& F;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_z;
        
        typedef typename Basis3D::FirstBasisType Basis_x;
        typedef typename Basis3D::SecondBasisType Basis_y;
        typedef typename Basis3D::ThirdBasisType Basis_z;

        Integral<Gauss, Basis_x, Basis_x>  integralf_x;
        Integral<Gauss, Basis_y, Basis_y> integralf_y;
        Integral<Gauss, Basis_z, Basis_z>  integralf_z;

    public:
        SeparableRHS3D(const Basis3D& _basis, const SeparableFunction3D<T>& _F,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y,
                       const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_z,
                       int order);

        T
        operator()(XType xtype_x, int j_x, long k_x,
                   XType xtype_y, int j_y, long k_y,
                   XType xtype_z, int j_z, long k_z) const;

        T
        operator()(const Index3D &index) const;
};

} // namespace lawa

#include <lawa/righthandsides/separablerhs3d.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H
