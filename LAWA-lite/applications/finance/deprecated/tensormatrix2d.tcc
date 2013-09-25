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


/*
class TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T,Basis>,
                     Compression, LeftPreconditioner, RightPreconditioner>
*/
template <typename T, typename Basis, typename CGMYOperator,
          typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::TensorMatrix2D
               (const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator> &_a_operator,
                const SpaceTimeInitialCondition1D<T, Basis> &_a_initcond,
                const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right, Compression &_c,
                T entrybound, int NumOfRows, int NumOfCols)
    : a_operator(_a_operator), a_initcond(_a_initcond), p_left(_p_left), p_right(_p_right), c(_c),
      c_t(a_operator.basis.first), pde_c_x(a_operator.basis.second), cgmy_c_x(a_operator.basis.second,a_operator.cgmy_x.cgmy.Y),
      data_d_t(a_operator.d_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_id_t(a_operator.id_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_id_x(a_operator.id_x, pde_c_x, entrybound, NumOfRows, NumOfCols),
      data_cgmy_x(a_operator.cgmy_x, cgmy_c_x, entrybound, NumOfRows, NumOfCols)

{
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::operator()(const Index2D &row_index,
                                                                                 const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;

    T prec = 1.;
    const_coeff_it it_P_left_end  = P_left_data.end();
    const_coeff_it it_row_index   = P_left_data.find(row_index);
    if (it_row_index != it_P_left_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp = p_left(row_index);
        P_left_data[row_index] = tmp;
        prec *= tmp;
    }
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }
    return prec * ( data_d_t(row_index.index1,col_index.index1) * data_id_x(row_index.index2,col_index.index2) +
                    data_id_t(row_index.index1,col_index.index1) * data_cgmy_x(row_index.index2,col_index.index2)
                  );
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>, Compression,
               LeftPreconditioner, RightPreconditioner>::operator()(const Index1D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }

    return prec * a_initcond(row_index,col_index);
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::left_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_left_end       = P_left_data.end();
    const_coeff_it it_index   = P_left_data.find(index);
    if (it_index != it_P_left_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_left(index);
        P_left_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::right_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end       = P_right_data.end();
    const_coeff_it it_index   = P_right_data.find(index);
    if (it_index != it_P_right_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_right(index);
        P_right_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
void
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::clear()
{
    data_d_t.clear();
    data_id_t.clear();
    data_cgmy_x.clear();
    data_id_x.clear();
}

}    //namespace lawa
