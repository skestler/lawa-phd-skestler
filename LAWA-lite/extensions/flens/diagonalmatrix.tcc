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

#include <cassert>

namespace flens {

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(const DenseVector<Array<T> > &diag)
    : _diag(diag)
{
}

template <typename T>
const DenseVector<Array<T> > &
DiagonalMatrix<T>::diag() const
{
    return _diag;
}

template <typename T>
int
DiagonalMatrix<T>::numRows() const
{
    return _diag.length();
}

template <typename T>
int
DiagonalMatrix<T>::numCols() const
{
    return _diag.length();
}

template <typename ALPHA,typename T,typename VX,typename BETA,typename VY>
void
mv(cxxblas::Transpose /*trans*/, ALPHA alpha, const DiagonalMatrix<T> &D,
   const VX &x, BETA beta, VY &y)
{
    assert(D.diag().length()==x.length());
    const DenseVector<Array<T> > &diag = D.diag();

    if (beta==0) {
        y.engine().resize(x.range()) || y.engine().fill();
    } else {
        y *= beta;
    }
    int cntr1=y.firstIndex(), cntr2=diag.firstIndex();
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i, ++cntr1, ++cntr2) {
        y(cntr1) += alpha * diag(cntr2)*x(i);
    }
}

template <typename T>
const DiagonalMatrix<T>
inv(const DiagonalMatrix<T> &D)
{
    const DenseVector<Array<T> > &diag = D.diag();
    DenseVector<Array<T> > v(diag.range());
    for (int i=diag.firstIndex(); i<=diag.lastIndex(); ++i) {
        v(i)=1.0/diag(i);
    }
    return DiagonalMatrix<T>(v);
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const DiagonalMatrix<T> &D)
{
    const DenseVector<Array<T> > & diag = D.diag();
    out << "diag( ";
    for (int i=diag.firstIndex(); i<diag.lastIndex(); ++i) {
        out << diag(i) << ", ";
    }
    out << diag(diag.lastIndex()) << " )";
    return out;
}


} // namespace flens

