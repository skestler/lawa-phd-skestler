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


#ifndef LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_SPARSEGRID2D_H
#define LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_SPARSEGRID2D_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/operators/operators.h>


namespace lawa {

template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
class SparseGrid2D{
        typedef IndexSet<Index2D>::const_iterator const_set_it;

    public:
        SparseGrid2D(const Basis2D &_basis, BilinearForm &_a, RHS &_rhsintegral, T _J, T _gamma);

        IndexSet<Index2D>
        getIndexSet();

        void
        solve_cg(T &energynorm=0);

        T
        evaluate(T x, T y);

    private:
        void
        setupIndexSet();

        const Basis2D        &basis;
        BilinearForm   &a;
        RHS            &rhs;
        DenseVector<Array<T> > u;
        T J;
        T gamma;
        IndexSet<Index2D>    Lambda;
};

}   //namespace lawa

#include <lawa/methods/sparsegrid/datastructures/sparsegrid2d.tcc>

#endif  // LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_SPARSEGRID2D_H
