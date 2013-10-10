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

#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS1D_WEIGHTEDSOBOLEVNORMPRECONDITIONER1D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS1D_WEIGHTEDSOBOLEVNORMPRECONDITIONER1D_H 1

#include <lawa/integrals/integrals.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>
#include <lawa/aux/compiletime_assert.h>

/*
 * Preconditioner for the space H^1_w := \{ v \in L_{1,loc}: \int v^2 w dx,
 *                                                           \int v_x^2 w dx < \infty \}
 */

namespace lawa {

template <typename T, typename Basis>
class WeightedSobolevNormPreconditioner1D
{
    ct_assert(IsPrimal<Basis>::value);

    public:
        WeightedSobolevNormPreconditioner1D(const Basis &basis, const Function<T> &weight,
                                            const int sobolev_order);

        T
        operator()(XType xtype1, int j1, long k1) const;

        T
        operator()(const Index1D &index) const;

    private:
        const int _sobolev_order;
        IntegralF<Gauss, Basis, Basis> _integral;
};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners1d/weightedsobolevnormpreconditioner1d.tcc>

#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS1D_WEIGHTEDSOBOLEVNORMPRECONDITIONER1D_H

