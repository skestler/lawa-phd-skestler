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

#ifndef LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMTENSORMATRIX2D_H
#define LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMTENSORMATRIX2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/uniform/algorithms/assembler1d.h>

namespace lawa {

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
         typename S2_x, typename S2_y>
class UniformTensorMatrix2D
{
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

    public:
        UniformTensorMatrix2D(const UniformBasis2D& basis, const S1_x &s1_x, const S1_y &s1_y,
                        const S2_x &s2_x, const S2_y &s2_y, const int _Jx, const int _Jy);

        int
        numRows() const;

        int
        numCols() const;

        IndexSet<Index2D>
        getIndexSet() const;

        DenseVectorT
        operator*(const DenseVectorT &v) const;

    private:
        void
        assembleMatrices();

        DenseMatrixT
        operator*(const DenseMatrixT &V) const;

        const UniformBasis2D &_basis;
        const S1_x  &_s1_x;
        const S1_y  &_s1_y;
        const S2_x  &_s2_x;
        const S2_y  &_s2_y;
        const int   _Jx, _Jy;

        Assembler1D<T,typename UniformBasis2D::FirstBasisType>  assembler_x;
        Assembler1D<T,typename UniformBasis2D::SecondBasisType> assembler_y;

        SparseMatrixT _M_s1_x, _M_s1_y, _M_s2_x, _M_s2_y;
};


} // namespace lawa

#include <lawa/methods/uniform/datastructures/uniformtensormatrix2d.tcc>

#endif // LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMTENSORMATRIX2D_H

