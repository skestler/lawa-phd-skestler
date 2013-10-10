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


#ifndef LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H
#define LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/uniform/algorithms/blockassembler1d.h>
#include <lawa/operators/operators.h>


/*
  std::vector<int*> sg_blocks:
     stores an array of "blocks" (see below). Each block consists of a level combination w.r.t
     to the sparse grid structure.



  int* pair:
     (pair[0], pair[1]) combination of levels in x and y direction w.r.t. to sparse grid
     structure,
     (pair[2], pair[3]) give the first and last index for the the current block
     when all sparse grid blocks are stored in one vector
     (pair[4], pair[5]) give the dimension of the block, i.e. pair[4] = |# nof x-indices|,
     pair[5] = |# nof y-indices|.
 */


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

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
class TensorSparseGrid2D {

    public:
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
        typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
        typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

        typedef std::map<std::pair<int,int>, int ,lt_int_vs_int >           LevelPairMap;

        TensorSparseGrid2D(const Basis2D &basis, const S1_x &s1_x, const S1_y &s1_y,
                           const S2_x &s2_x, const S2_y &s2_y, int I, T gamma);

        int
        getDimension() const;

        int
        numCols() const;

        int
        numRows() const;

        IndexSet<Index2D>
        getIndexSet() const;

        void
        toCoefficients(const DenseVectorT &vec,
                       Coefficients<Lexicographical,T,Index2D> &sparsegridcoefficients);

        T
        evaluate(const DenseVectorT &u, T x, T y, int deriv_x, int deriv_y) const;

        DiagonalMatrixT
        assembleDiagonalMatrixPreconditioner();

        template <typename RHSIntegral_x, typename RHSIntegral_y>
        DenseVectorT
        assembleRHS(const RHSIntegral_x &rhs_x, const RHSIntegral_y &rhs_y);

        DenseVectorT
        operator*(const DenseVectorT &v) const;

    private:
        DenseMatrixT
        block_multiplication(int row_block, int col_block, const DenseMatrixT &Xj) const;

        void
        assembleMatrices();

        const Basis2D       &_basis;
        const S1_x          &_s1_x;
        const S1_y          &_s1_y;
        const S2_x          &_s2_x;
        const S2_y          &_s2_y;
        int                 _j0_x, _j0_y;
        int                 _I;
        T                   _gamma;

        BlockAssembler1D<T,typename Basis2D::FirstBasisType> _blockassembler1d;
        int                                                  _dim;
        std::vector<int*>                                    _sg_blocks;
        LevelPairMap                                         _levelpair_map;
        std::vector<SparseMatrixT>                           _matrixblocks_s1_x, _matrixblocks_s1_y,
                                                             _matrixblocks_s2_x, _matrixblocks_s2_y;
};

}   //namespace lawa

#include <lawa/methods/sparsegrid/datastructures/tensorsparsegrid2d.tcc>

#endif  // LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H
