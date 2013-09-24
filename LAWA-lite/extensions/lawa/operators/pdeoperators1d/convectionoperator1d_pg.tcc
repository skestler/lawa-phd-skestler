/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

namespace lawa {


template <typename T, typename TrialBasis, typename TestBasis>
ConvectionOperator1D_PG<T, TrialBasis, TestBasis>::ConvectionOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis)
  : trialbasis(_trialbasis), testbasis(_testbasis), integral(_testbasis, _trialbasis)
{
}

template <typename T, typename TrialBasis, typename TestBasis>
T
ConvectionOperator1D_PG<T, TrialBasis, TestBasis>::operator()(XType xtype1, int j1, long k1,
                                                              XType xtype2, int j2, long k2) const
{   
    // v * u_x
    return integral(j1, k1, xtype1, 0, j2, k2, xtype2, 1);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
ConvectionOperator1D_PG<T, TrialBasis, TestBasis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return ConvectionOperator1D_PG<T, TrialBasis, TestBasis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                                      col_index.xtype, col_index.j, col_index.k);
}



}    //namespace lawa
