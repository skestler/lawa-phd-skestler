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

namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::MapMatrix(const BilinearForm &_a,
                                                                      const Preconditioner &_p,
                                                                      Compression &_compression)
:   data(6291469), a(_a), p(_p), compression(_compression), P_data()
{
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(const Index &row_index,
                                                                       const Index &col_index)
{
    Entry<Index> entry(row_index,col_index);

    typedef typename EntryMap::const_iterator const_map_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    const_map_it it_end = data.end();
    const_map_it it_entry = data.find(entry);

    if (it_entry != it_end) {
        return (*it_entry).second;
    }
    else {
        T prec = 1.;
        if (!flens::IsSame<NoPreconditioner<T, Index>, Preconditioner>::value) {
            const_coeff_it it_P_end       = P_data.end();
            const_coeff_it it_row_index = P_data.find(row_index);
            if (it_row_index != it_P_end) {
                prec *= (*it_row_index).second;
            }
            else {
                T tmp = p(row_index);
                P_data[row_index] = tmp;
                prec *= tmp;
            }
            it_P_end       = P_data.end();
            const_coeff_it it_col_index   = P_data.find(col_index);
            if (it_col_index != it_P_end) {
                prec *= (*it_col_index).second;
            }
            else {
                T tmp = p(col_index);
                P_data[col_index] = tmp;
                prec *= tmp;
            }
        }
        T val = prec * a(row_index,col_index);

        if (fabs(val) > 0.) data.insert(val_type(entry,val));
        return val;
    }

}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
Coefficients<Lexicographical,T,Index>
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::mv(const IndexSet<Index> &LambdaRow,
                                                    const Coefficients<Lexicographical,T,Index> &x)
{
    return lawa::mv(LambdaRow,*this,x);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix
                                (const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                                 flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens,
                                 int /*J*/)
{
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix
                                (const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                                 flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens,
                                 T /*tol*/)
{
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename Index, typename BilinearForm, typename Compression,
          typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::
apply(const Coefficients<Lexicographical,T,Index> &u, T tol, const IndexSet<Index> Lambda,
      Coefficients<Lexicographical,T,Index> &r)
{
    std::cerr << "Apply not implemented for MapMatrix." << std::endl;
    exit(1);
}
/*
template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(T t, const  Index &row_index, const Index &col_index)
{
    Entry<Index> entry(row_index,col_index);

    typedef typename EntryMap::const_iterator const_map_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    const_map_it it_end = data.end();
    const_map_it it_entry = data.find(entry);

    if (it_entry != it_end) {
         return (*it_entry).second;
    }
    else {
          T prec = 1.;
           const_coeff_it it_P_end       = P_data.end();
           const_coeff_it it_row_index = P_data.find(row_index);
           if (it_row_index != it_P_end) {
               prec *= (*it_row_index).second;
           }
           else {
               T tmp = p(row_index);
               P_data[row_index] = tmp;
               prec *= tmp;
           }
           it_P_end       = P_data.end();
           const_coeff_it it_col_index   = P_data.find(col_index);
           if (it_col_index != it_P_end) {
               prec *= (*it_col_index).second;
           }
           else {
               T tmp = p(col_index);
               P_data[col_index] = tmp;
               prec *= tmp;
           }
           T val = prec * a(row_index,col_index);
           if (fabs(val) > 0) data.insert(val_type(entry,val));
           return val;
    }
}
*/

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::clear()
{
    data.clear();
}

} // namespace lawa

