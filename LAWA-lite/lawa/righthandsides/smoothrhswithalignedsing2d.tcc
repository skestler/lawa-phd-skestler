/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

template<typename T, typename Basis2D, QuadratureType Quad>
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::SmoothRHSWithAlignedSing2D
                                            (const Basis2D &_basis, const Function2D<T> &_F,
                                             int order, unsigned short _derivx,
                                             unsigned short _derivy)
   : basis(_basis), integral2d(_F, basis.first, basis.second), derivx(_derivx), derivy(_derivy)
{
    integral2d.quadrature.setOrder(order);
}

template<typename T, typename Basis2D, QuadratureType Quad>
T
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::operator()(XType ex, int jx, long kx,
                                                       XType ey, int jy, long ky) const
{
    return integral2d(jx, kx, ex, derivx, jy, ky, ey, derivy);
}

template<typename T, typename Basis2D, QuadratureType Quad>
T
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}    //namespace lawa
