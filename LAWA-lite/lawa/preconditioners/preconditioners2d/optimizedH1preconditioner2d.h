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

#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS2D_OPTIMIZEDH1PRECONDITIONER2D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS2D_OPTIMIZEDH1PRECONDITIONER2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>
#include <lawa/aux/compiletime_assert.h>

namespace lawa {

template <typename T, typename Basis2D>
class OptimizedH1Preconditioner2D
{
    public:
        OptimizedH1Preconditioner2D(const Basis2D &_basis, T _a_x=1., T _a_y=1., T _c=1.);

        void
        setParameters(T _a_x, T _a_y, T _c);

        void
        setThetaTimeStepParameters(T theta, T timestep);

        T
        operator()(XType xtype_x, int j_x, long k_x,
                   XType xtype_y, int j_y, long k_y) const;

        T
        operator()(const Index2D &index) const;

        T
        operator[](const Index2D &index) const;

    private:
        typedef typename Basis2D::FirstBasisType Basis_x;
        typedef typename Basis2D::SecondBasisType Basis_y;

        const Basis_x &basis_x;
        const Basis_y &basis_y;

        T a_x, a_y, c;
        T factor_a_x, factor_a_y, factor_c;
};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners2d/optimizedH1preconditioner2d.tcc>

#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS2D_OPTIMIZEDH1PRECONDITIONER2D_H

