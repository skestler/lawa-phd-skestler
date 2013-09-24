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
 
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename Basis, typename Compression, typename Preconditioner>
AdaptiveHelmholtzOperator3D<T, Basis, HelmholtzOperator3D<T, Basis>, Compression, Preconditioner>::
AdaptiveHelmholtzOperator3D(const Basis3D &_basis, T _c, const Preconditioner &_Prec, 
                            T entrybound=0., int NumOfRows=4096, int NumOfCols=2048)
    : basis(_basis), c(_c), Prec(_Prec),
      compression_1d_x(_basis.first), compression_1d_y(_basis.second), compression_1d_z(_basis.third),
      op_identity_x(_basis.first), op_identity_y(_basis.second), op_identity_z(_basis.third),
      op_laplace_x(_basis.first), op_laplace_y(_basis.second), op_laplace_z(_basis.third),
      entrybound(_entrybound), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
      data_identity_x(op_identity_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
      data_identity_y(op_identity_y, compression_1d_y, entrybound, NumOfRows, NumOfCols),
      data_identity_z(op_identity_z, compression_1d_z, entrybound, NumOfRows, NumOfCols),
      data_laplace_x(op_laplace_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
      data_laplace_y(op_laplace_y, compression_1d_y, entrybound, NumOfRows, NumOfCols),
      data_laplace_z(op_laplace_z, compression_1d_z, entrybound, NumOfRows, NumOfCols),
      P_data()
{
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
T
AdaptiveHelmholtzOperator3D<T, Basis, HelmholtzOperator3D<T, Basis>, Compression, Preconditioner>::
operator()(const Index3D &row_index, const Index3D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, Prec>::value) {
        // Row Precondioning:
        const_coeff_it it_row_index   = P_data.find(row_index);
        //  Entry has already been computed:
        if (it_row_index != P_data.end()) {
            prec *= (*it_row_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = Prec(row_index);
            P_data[row_index] = tmp;
            prec *= tmp;
        }
    }
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, Prec>::value) {
        // Col precondioning:
        const_coeff_it it_col_index   = P_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = Prec(col_index);
            P_data[col_index] = tmp;
            prec *= tmp;
        }
    }
    

    T id_x = data_identity_x(row_index.index1,col_index.index1);
    T id_y = data_identity_y(row_index.index2,col_index.index2);
    T id_z = data_identity_z(row_index.index3,col_index.index3);
    
    T dd_x = data_laplace_x(row_index.index1,col_index.index1);
    T dd_y = data_laplace_y(row_index.index2,col_index.index2);
    T dd_z = data_laplace_z(row_index.index3,col_index.index3);

    return prec *( dd_x*id_y*id_z + id_x*dd_y*id_z + id_x*id_y*dd_z + id_x*id_y*id_z );
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
T
AdaptiveHelmholtzOperator3D<T, Basis, HelmholtzOperator3D<T, Basis>, Compression, Preconditioner>::prec(const Index3D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff_it;
    T prec = 1.;
   
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, Prec>::value) {
        const_coeff_it it_index   = P_data.find(index);
        //  Entry has already been computed:
        if (it_index != P_data.end()) {
            prec *= (*it_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = Prec(index);
            P_data[index] = tmp;
            prec *= tmp;
        }
    }

    return prec;
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
void
AdaptiveHelmholtzOperator3D<T, Basis, HelmholtzOperator3D<T, Basis>, Compression, Preconditioner>::clear()
{
    data_identity_x.clear();
    data_identity_y.clear();
    data_identity_z.clear();
    data_laplace_x.clear();
    data_laplace_y.clear();
    data_laplace_y.clear();
}

}    //namespace lawa

