/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef EXTENSIONS_FLENS_DIAGONALMATRIX_H
#define EXTENSIONS_FLENS_DIAGONALMATRIX_H 1

#include <iostream>
#include <lawa/flensforlawa.h>

namespace flens {

template <typename T>
struct DiagonalMatrix;

template <typename T>
struct TypeInfo<DiagonalMatrix<T> >
{
    typedef DiagonalMatrix<T> Impl;
    typedef T                 ElementType;
};

template <typename T>
class DiagonalMatrix : 
    public SymmetricMatrix<DiagonalMatrix<T> >
{
    public:
        typedef T                 ElementType;

        DiagonalMatrix(const DenseVector<Array<T> > &diag);

        const DenseVector<Array<T> > &
        diag() const;

        int
        numRows() const;
    
        int
        numCols() const;
    
//    private:
        DenseVector<Array<T> > _diag;
};

//------------------------------------------------------------------------------

template <typename ALPHA,typename T,typename VX,typename BETA,typename VY>
    void
    mv(cxxblas::Transpose trans, ALPHA alpha, const DiagonalMatrix<T> &D,
    const VX &x, BETA beta, VY &y);


template <typename T>
    const DiagonalMatrix<T>
    inv(const DiagonalMatrix<T> &D);

template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const DiagonalMatrix<T> &D);

} // namespace flens

#include <extensions/flens/diagonalmatrix.tcc>

#endif // EXTENSIONS_FLENS_DIAGONALMATRIX_H 

