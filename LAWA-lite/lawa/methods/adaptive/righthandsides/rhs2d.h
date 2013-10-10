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
#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS2D_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS2D_H 1

#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
struct RHS2D
{
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;

    RHS2D(const RHSINTEGRAL &rhsintegral, Preconditioner &P);

    bool
    readIndexSets(const char *filename);

    T
    operator()(const Index2D &lambda);

    Coefficients<Lexicographical,T,Index2D>
    operator()(const IndexSet<Index2D> &Lambda);

    Coefficients<Lexicographical,T,Index2D>
    operator()(T tol);

    IndexSet<Index2D>
    getFullIndexSet();

    void
    initializePropagation(const Coefficients<Lexicographical,T,Index2D> &f)
    {
        std::cerr << "RHS2D::initializePropagation is a dummy routine and should not be called" << std::endl;
        exit(1);
    }

    const RHSINTEGRAL                       &rhsintegral;
    Preconditioner                          &P;
    Coefficients<Lexicographical,T,Index2D> rhs_data;
    std::list<IndexSet<Index2D> >           rhs_indexsets;
    T                                       norm_estimate;
    T                                       current_tol;
    long double                             current_ell2norm;
    bool                                    new_values;
};

}    //namespace lawa

#include <lawa/methods/adaptive/righthandsides/rhs2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_RHS1D_H
