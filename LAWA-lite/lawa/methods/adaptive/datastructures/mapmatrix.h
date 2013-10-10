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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPMATRIX_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPMATRIX_H 1

#include <utility>
#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#else
    #include <ext/hash_map>
#endif
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/aux/timer.h>

namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
struct MapMatrix
{
    //typedef typename std::map<Entry<Index>,T,lt<Lexicographical,Index > > EntryMap;
#ifdef TRONE
    typedef typename std::tr1::unordered_map<Entry<Index>, T, entry_hashfunction<Index>,
                                                             entry_eqfunction<Index> > EntryMap;
#elif BOOST
    typedef typename boost::unordered_map<Entry<Index>, T, entry_hashfunction<Index>,
                                                                 entry_eqfunction<Index> > EntryMap;
#else
    typedef typename __gnu_cxx::hash_map<Entry<Index>, T, entry_hashfunction<Index>,
                                                             entry_eqfunction<Index> > EntryMap;
#endif
    typedef typename EntryMap::value_type val_type;
    EntryMap data;

    const BilinearForm &a;
    const Preconditioner &p;
    Compression &compression;
    Coefficients<Lexicographical,T,Index> P_data;

    MapMatrix(const BilinearForm &a, const Preconditioner &p, Compression &_compression);

    T
    operator()(const Index &row_index, const Index &col_index);

    Coefficients<Lexicographical,T,Index>
    mv(const IndexSet<Index> &LambdaRow,
       const Coefficients<Lexicographical,T,Index> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens, int J=-1);

    void
    toFlensSparseMatrix(const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens, T tol);

    void
    apply(const Coefficients<Lexicographical,T,Index> &u, T tol, const IndexSet<Index> Lambda,
          Coefficients<Lexicographical,T,Index> &r);

    //T
    //operator()(T t, const  Index &row_index, const Index &col_index);

    void
    clear();
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/mapmatrix.tcc>

#endif //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPMATRIX_H

