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


#ifndef  APPLICATIONS_FINANCE_TENSORMATRIX2D_H
#define  APPLICATIONS_FINANCE_TENSORMATRIX2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/constructions/tensorbasis.h>
#include <lawa/operators/operators.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/compressions/compressions.h>
#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>


template <typename T, typename Basis, typename OperatorBilinearForm, typename OperatorInitCond,
          typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
struct TensorMatrix2D {
};

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
struct TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T,Basis>,
                     Compression, LeftPreconditioner, RightPreconditioner>
{
    typedef CompressionPDE1D<T, Index1D, typename Basis::FirstBasisType>               Compression_t;
    typedef CompressionPDE1D<T, Index1D, typename Basis::SecondBasisType>            PDECompression_x;
    typedef CompressionCGMYOperator1D<T, typename Basis::SecondBasisType> CGMYCompression_x;

    typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;

    typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Convection_t, Compression_t, NoPreconditioner1D> DataConvection_t;
    typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Reaction_t,  Compression_t, NoPreconditioner1D>  DataReaction_t;
    typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Reaction_x,  PDECompression_x, NoPreconditioner1D>  DataReaction_x;
    typedef MapMatrixWithZeros<T, Index1D, CGMYOperator, CGMYCompression_x, NoPreconditioner1D>  DataCGMY_x;

    const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>     &a_operator;
    const SpaceTimeInitialCondition1D<T, Basis> &a_initcond;
    const LeftPreconditioner  &p_left;
    const RightPreconditioner &p_right;
    Compression &c;
    Coefficients<Lexicographical,T,Index2D> P_left_data;
    Coefficients<Lexicographical,T,Index2D> P_right_data;

    Compression_t     c_t;
    PDECompression_x  pde_c_x;
    CGMYCompression_x cgmy_c_x;

    DataConvection_t data_d_t;
    DataReaction_t   data_id_t;
    DataCGMY_x       data_cgmy_x;
    DataReaction_x   data_id_x;

    TensorMatrix2D(const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator> &_a_operator,
                   const SpaceTimeInitialCondition1D<T,Basis> &_a_initcond,
                   const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right,
                   Compression &c, T _entrybound=0., int NumOfRows=4096, int NumOfCols=2048);

    //call of p_left * a_operator * p_right
    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    //call of a_initcond * p_right
    T
    operator()(const Index1D &row_index, const Index2D &col_index);

    T
    left_prec(const Index2D &index);

    T
    right_prec(const Index2D &index);

    void
    clear();

};

#include <applications/finance/tensormatrix2d.tcc>

#endif // APPLICATIONS_FINANCE_TENSORMATRIX2D_H
