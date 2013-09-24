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

#include <cmath>

namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::MapMatrixWithZeros
                                                                    (const BilinearForm &_a,
                                                                     const Preconditioner &_p,
                                                                     Compression &_compression,
                                                                     T _entrybound,
                                                                     int _NumOfRows, int _NumOfCols)
:  a(_a), p(_p), compression(_compression), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
   ConsecutiveIndices(), Zeros( (NumOfRows*NumOfCols) >> 5), warning_overflow(false),
   entrybound(_entrybound)
{
    PrecValues.engine().resize(int(NumOfRows));
    Zeros.assign((NumOfRows*NumOfCols) >> 5, (long) 0);
    //NonZeros.resize(3145739);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
T
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::operator()
                                                                     (const Index &_row_index,
                                                                      const Index &_col_index)
{
    Index temp_row_index = _row_index;
    Index temp_col_index = _col_index;
    typedef typename IndexSet<Index>::iterator set_it;
    int size = (int)ConsecutiveIndices.size();
    std::pair<set_it,bool> row = ConsecutiveIndices.insert(temp_row_index);
    std::pair<set_it,bool> col = ConsecutiveIndices.insert(temp_col_index);

    set_it row_index = row.first;
    if (row.second) {
        (*row_index).linearindex = size;
        ++size;
    }
    set_it col_index = col.first;
    if (col.second) {
        (*col_index).linearindex = size;
    }

    unsigned long value, block;
    unsigned int block_num, block_pos;

    if (((*row_index).linearindex >= NumOfRows) || ((*col_index).linearindex >= NumOfCols)) {
        value = 3;
    }
    else {
        block_num = (   (*col_index).linearindex*NumOfRows + (*row_index).linearindex ) >> 5;
        block_pos = ( ( (*col_index).linearindex*NumOfRows + (*row_index).linearindex ) & 31 ) * 2;

        block = Zeros[block_num];
        //long value = ( (((long) 3) << (62-block_pos)*2) & (block) ) >> (62-block_pos);
        value = ( (((long) 3) << block_pos) & (block) ) >> block_pos;
    }


    if (value == 1) {
        return 0.;
    }
    else if (value == 2) {
        T tmp = NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)];
        return tmp;
    }
    else if (value == 0) {     //value == 0
        T prec = 1.;
        if (!flens::IsSame<NoPreconditioner<T, Index>, Preconditioner>::value) {
            if (fabs(PrecValues((*row_index).linearindex+1)) > 0) {
                prec *= PrecValues((*row_index).linearindex+1);
            }
            else {
                T tmp = p(*row_index);
                prec *= tmp;
                PrecValues((*row_index).linearindex+1) = tmp;
            }
            if (fabs(PrecValues((*col_index).linearindex+1)) > 0) {
                prec *= PrecValues((*col_index).linearindex+1);

            }
            else {
                T tmp = p(*col_index);
                prec *= tmp;
                PrecValues((*col_index).linearindex+1) = tmp;
            }
        }
        T val = 0.;
        val = prec * a(*row_index,*col_index);
        if (fabs(val)>entrybound) {
            NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)] = val;
            Zeros[block_num] = (((long) 2) << block_pos) | (block) ;
            return val;
        }
        else {
            Zeros[block_num] = (((long) 1) << block_pos) | (block) ;
            return 0.;
        }
    }
    else {
        if (!warning_overflow) {
            std::cout << "WARNING: Not enough storage for zeros!!" << std::endl;
            warning_overflow = true;
        }
        T prec = 1.;
        if (   ((*row_index).linearindex < NumOfRows)
            && (fabs(PrecValues((*row_index).linearindex+1)) > 0)) {
            prec *= PrecValues((*row_index).linearindex+1);
        }
        else {
            prec *= p(*row_index);
        }
        if (   ((*col_index).linearindex < NumOfRows)
            && (fabs(PrecValues((*col_index).linearindex+1)) > 0)) {
            prec *= PrecValues((*col_index).linearindex+1);
        }
        else {
            prec *= p(*col_index);
        }
        T val = 0.;
        val = prec * a(*row_index,*col_index);
        if (fabs(val)>entrybound) return val;
        else                       return 0.;
    }
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
Coefficients<Lexicographical,T,Index>
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::mv
                                                    (const IndexSet<Index> &LambdaRow,
                                                     const Coefficients<Lexicographical,T,Index> &x)
{
    return lawa::mv(LambdaRow,*this,x);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
void
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix
                                (const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                                 flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens,
                                 int /*J*/)
{
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
void
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix
                                (const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                                 flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens,
                                 T /*tol*/)
{
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
void
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::clear()
{
    NonZeros.clear();
    Zeros.clear();
    ConsecutiveIndices.clear();
    PrecValues.engine().resize(0);

}

}    //namespace lawa

