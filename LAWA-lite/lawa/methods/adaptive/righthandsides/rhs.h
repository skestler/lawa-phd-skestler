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

#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS_H 1

#include <lawa/methods/adaptive/algorithms/thresh.h>
#include <lawa/methods/adaptive/datastructures/adaptiverhs.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
class RHS : public AdaptiveRhs<T, Index>
{
    public:
        const RHSINTEGRAL &rhsintegral;
        Preconditioner &P;
        Coefficients<Lexicographical,T,Index> rhs_data;
        Coefficients<AbsoluteValue,T,Index>   rhs_abs_data;

    //public:
        RHS(const RHSINTEGRAL &rhsintegral, Preconditioner &P);

        RHS(const RHSINTEGRAL &rhsintegral, Preconditioner &P,
            const Coefficients<Lexicographical,T,Index> &_rhs_data);

        T
        operator()(const Index &lambda);

        Coefficients<Lexicographical,T,Index>
        operator()(const IndexSet<Index> &Lambda);

        Coefficients<Lexicographical,T,Index>
        operator()(T tol);

        T
        operator()(T t, const Index &lambda);

        Coefficients<Lexicographical,T,Index>
        operator()(T t, const IndexSet<Index> &Lambda);

};

}    //namespace lawa

#include <lawa/methods/adaptive/righthandsides/rhs.tcc>

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS_H

