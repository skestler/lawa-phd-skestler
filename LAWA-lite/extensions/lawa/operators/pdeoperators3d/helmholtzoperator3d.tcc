/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Kristina Steih, Mario Rometsch, Alexander Stippler.

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

namespace lawa {

template <typename T, typename Basis3D>
HelmholtzOperator3D<T, Basis3D>::HelmholtzOperator3D(const Basis3D & _basis, const T _c)
    : basis(_basis), c(_c),
      integral_x(_basis.first, _basis.first), integral_y(_basis.second, _basis.second),
      integral_z(_basis.third, _basis.third)
{
}

template <typename T, typename Basis3D>
T
HelmholtzOperator3D<T, Basis3D>::operator()(XType row_xtype_x, int j1_x, long k1_x,
                                            XType row_xtype_y, int j1_y, long k1_y,
                                            XType row_xtype_z, int j1_z, long k1_z,
                                            XType col_xtype_x, int j2_x, long k2_x,
                                            XType col_xtype_y, int j2_y, long k2_y,
                                            XType col_xtype_z, int j2_z, long k2_z) const
{
    T val_x =    integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0);
    T dd_val_x = integral_x(j1_x, k1_x, row_xtype_x, 1, j2_x, k2_x, col_xtype_x, 1);
    T val_y =    integral_y(j1_y, k1_y, row_xtype_y, 0, j2_y, k2_y, col_xtype_y, 0);
    T dd_val_y = integral_y(j1_y, k1_y, row_xtype_y, 1, j2_y, k2_y, col_xtype_y, 1);
    T val_z =    integral_z(j1_z, k1_z, row_xtype_z, 0, j2_z, k2_z, col_xtype_z, 0);
    T dd_val_z = integral_z(j1_z, k1_z, row_xtype_z, 1, j2_z, k2_z, col_xtype_z, 1);

    // (v1_x * u1_x)*(v2 * u2)*(v3 * u3) + (v1 * u1)*(v2_y * u2_y)*(v3 * u3) + (v1 * u1)*(v2 * u2)*(v3_x * u3_x)
    // + c * (v1 * u1)*(v2 * u2)*(v3 * u3)
    return dd_val_x * val_y * val_z + val_x * dd_val_y * val_z + val_x * val_y * dd_val_z + c * val_x * val_y * val_z;
}

template <typename T, typename Basis3D>
T
HelmholtzOperator3D<T, Basis3D>::operator()(const Index3D &row_index, const Index3D &col_index) const
{
    return HelmholtzOperator3D<T, Basis3D>::operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                                                       row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                                                       row_index.index3.xtype, row_index.index3.j, row_index.index3.k,
                                                       col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                                                       col_index.index2.xtype, col_index.index2.j, col_index.index2.k,
                                                       col_index.index3.xtype, col_index.index3.j, col_index.index3.k);
}

}    //namespace lawa

