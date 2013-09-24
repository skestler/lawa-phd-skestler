/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY1D_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
class SYM_APPLY_1D {

    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_coeff_abs_it;
    typedef typename IndexSet<Index>::const_iterator const_set_it;

    public:
        typedef MA MAType;

        const Parameters &parameters;
        const Basis1D &basis;
        MA &A;

        SYM_APPLY_1D(const Parameters &_parameters, const Basis1D &_basis, MA &_A);

        Coefficients<Lexicographical,T,Index>
        operator()(const Coefficients<Lexicographical,T,Index> &v, int k);

        Coefficients<Lexicographical,T,Index>
        operator()(const Coefficients<Lexicographical,T,Index> &v, int k, int J);

        Coefficients<Lexicographical,T,Index>
        operator()(const Coefficients<Lexicographical,T,Index> &v, T eps);

    private:

        int
        findK(const Coefficients<AbsoluteValue,T,Index> &v, T eps);

};

}   //namespace lawa

#include <lawa/methods/adaptive/algorithms/apply1d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_APPLY1D_H

