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


#ifndef LAWA_METHODS_UNIFORM_ALGORITHMS_BLOCKASSEMBLER1D_H
#define LAWA_METHODS_UNIFORM_ALGORITHMS_BLOCKASSEMBLER1D_H 1

#include <extensions/extensions.h>

namespace lawa{

template<typename T, typename Basis>
class BlockAssembler1D
{
    private:
        const Basis& basis;

    public:
        BlockAssembler1D(const Basis& _basis);

        template <typename BilinearForm>
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
        assembleStiffnessMatrixBlock(BilinearForm& a, int i1, int i2, T tol = 10e-15);

        template <typename BilinearForm>
        flens::DenseVector<flens::Array<T> >
        assembleStiffnessMatrixBlockDiagonal(BilinearForm& a, int i, T tol = 10e-15);

        template <typename RHSIntegral>
        flens::DenseVector<flens::Array<T> >
        assembleRHSBlock(RHSIntegral& rhs, int i);


};

} // namespace lawa

#include <lawa/methods/uniform/algorithms/blockassembler1d.tcc>

#endif // LAWA_METHODS_UNIFORM_ALGORITHMS_BLOCKASSEMBLER1D_H

