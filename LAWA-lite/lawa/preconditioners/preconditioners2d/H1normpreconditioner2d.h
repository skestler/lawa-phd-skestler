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

#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS2D_H1NORMPRECONDITIONER2D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS2D_H1NORMPRECONDITIONER2D_H 1

#include <lawa/integrals/integrals.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>
#include <lawa/aux/compiletime_assert.h>

namespace lawa {

template <typename T, typename Basis2D>
class H1NormPreconditioner2D
{
    public:
        H1NormPreconditioner2D(const Basis2D &_basis);

        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(const Index2D &index) const;

    private:
        typedef typename Basis2D::FirstBasisType Basis_x;
        typedef typename Basis2D::SecondBasisType Basis_y;

        Integral<Gauss, Basis_x, Basis_x>   _integral_x;
        Integral<Gauss, Basis_y, Basis_y>   _integral_y;

        T refval_id_bspline;
        T refval_dd_bspline;
        T refval_id_wavelet;
        T refval_dd_wavelet;
};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners2d/H1normpreconditioner2d.tcc>

#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS2D_H1NORMPRECONDITIONER2D_H

