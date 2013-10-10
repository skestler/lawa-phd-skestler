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


#ifndef LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H
#define LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/integrals/integral.h>

namespace lawa {

/* Helmholtz Operator 3D
 * 
 *      a(v,u) =      Integral(v1_x * u1_x) * Integral(v2 * u2)     * Integral(v3 * u3) 
 *              +     Integral(v1 * u1)     * Integral(v2_y * u2_y) * Integral(v3 * u3) 
 *              +     Integral(v1 * u1)     * Integral(v2 * u2)     * Integral(v3_x * u3_x)
 *              + c * Integral(v1 * u1)     * Integral(v2 * u2)     * Integral(v3 * u3)
 *
 */
template <typename T, typename Basis3D>
class HelmholtzOperator3D{

    private:

        typedef typename Basis3D::FirstBasisType  Basis_x;
        typedef typename Basis3D::SecondBasisType Basis_y;
        typedef typename Basis3D::ThirdBasisType  Basis_z;
        
        Integral<Gauss, Basis_x, Basis_x> integral_x;
        Integral<Gauss, Basis_y, Basis_y> integral_y;
        Integral<Gauss, Basis_z, Basis_z> integral_z;

    public:
        
        const Basis3D &basis;
        const T c;
        
        HelmholtzOperator3D(const Basis3D& _basis, const T _c);

        T
        operator()(XType row_xtype_x, int j1_x, long k1_x,
                   XType row_xtype_y, int j1_y, long k1_y,
                   XType row_xtype_z, int j1_z, long k1_z,
                   XType col_xtype_x, int j2_x, long k2_x,
                   XType col_xtpye_y, int j2_y, long k2_y,
                   XType col_xtpye_z, int j2_z, long k2_z) const;

        T
        operator()(const Index3D &row_index, const Index3D &col_index) const;
};

}    // namespace lawa

#include <lawa/operators/pdeoperators3d/helmholtzoperator3d.tcc>

#endif //  LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H

