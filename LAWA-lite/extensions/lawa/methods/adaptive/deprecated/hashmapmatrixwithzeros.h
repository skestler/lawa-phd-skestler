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


#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HASHMAPMATRIXWITHZEROS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HASHMAPMATRIXWITHZEROS_H 1

#include <utility>
#include <ext/hash_map>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/aux/timer.h>
#include <lawa/preconditioners/preconditioners.h>

#define ROW_SIZE_2D 4*8192
#define COL_SIZE_2D 2*4*2048

namespace lawa {

struct lt_int_vs_int
{
    inline
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) const
    {
        if (left.first != right.first) return left.first < right.first;
        else                           return left.second < right.second;
    }
};

struct hash_pair_of_int {
    inline
    size_t operator()(const std::pair<int,int>& p) const {
        return ( (p.first+p.second)*(p.first+p.second+1)/2 + p.second ) %  9369319;
    }
};

struct equal_pair_of_int {
    inline
    bool operator()(const std::pair<int,int>& p_left, const std::pair<int,int>& p_right) const {
        if (p_left.first != p_right.first) return false;
        else                               return (p_left.second == p_right.second);
    }
};

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
struct MapMatrixWithZeros
{
    typedef typename std::map<std::pair<int,int>,T,lt_int_vs_int > EntryMap;
    //typedef typename __gnu_cxx::hash_map<std::pair<int,int>, T, hash_pair_of_int,
    // equal_pair_of_int> EntryMap;
    typedef typename EntryMap::value_type val_type;

    EntryMap NonZeros;

    const BilinearForm &a;
    const Preconditioner &p;
    Compression &compression;

    unsigned int NumOfRows, NumOfCols;
    IndexSet<Index> ConsecutiveIndices;
    flens::DenseVector<Array<T> > PrecValues;
    std::vector<unsigned long> Zeros;
    bool warning_overflow;
    T entrybound;


    MapMatrixWithZeros(const BilinearForm &a, const Preconditioner &p, Compression &_compression,
                       T _entrybound=0., int NumOfRow=ROW_SIZE_2D, int NumOfCols=COL_SIZE_2D);

    T
    operator()(const Index &row_index, const Index &col_index);

    Coefficients<Lexicographical,T,Index>
    mv(const IndexSet<Index> &LambdaRow,
       const Coefficients<Lexicographical,T,Index> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens, int J);

    void
    toFlensSparseMatrix(const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens, T tol);

    void
    clear();
};


}    //namespace lawa

#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HASHMAPMATRIXWITHZEROS_H

