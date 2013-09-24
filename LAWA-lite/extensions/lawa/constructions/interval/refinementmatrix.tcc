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
#include <lawa/math/math.h>

namespace flens {

using namespace lawa;

template <typename T, Construction Cons>
RefinementMatrix<T,Interval,Cons>::RefinementMatrix()
{
}

template <typename T, Construction Cons>
RefinementMatrix<T,Interval,Cons>::RefinementMatrix(
                                int nLeft, int nRight,
                                const GeMatrix<FullStorage<T, ColMajor> > &A,
                                int _min_j0, int cons_j)
    : left(nLeft, (nLeft>0) ? A.firstCol() : 1),
      right(nRight, (nRight>0) ? A.lastCol()-nRight+1 : A.lastCol()+2),
      lengths(_(-nRight, nLeft)),
      min_j0(_min_j0), _cons_j(cons_j), _j(cons_j),
      _firstRow(A.firstRow()), _firstCol(A.firstCol()),
      _lastRow(A.lastRow()), _lastCol(A.lastCol()),
      _additionalRows(0), _additionalCols(0)
{
    assert(nLeft>=0);
    assert(nRight>=0);
    assert(_cons_j>=min_j0);
    
    assert(_firstCol>=1);

    _extractMasks(A);

    for (int i=left.firstIndex(); i<=left.lastIndex(); ++i) {
        lengths(1+i-left.firstIndex()) = left(i).length()+(A.firstRow()-1);
    }

    lengths(0) = leftband.firstIndex() - 1;//-A.firstRow();

    for (int i=right.firstIndex(); i<=right.lastIndex(); ++i) {
        lengths(-nRight+i-right.firstIndex()) = right(i).length()+(A.firstRow()-1);
    }
}

template <typename T, Construction Cons>
const typename DenseVector<Array<T> >::ConstView
RefinementMatrix<T,Interval,Cons>::operator()(int j, const Underscore<int> &/*u*/,
                                              int col) const
{
    assert(j>=min_j0);

    int additionalCols = 0, additionalRows = 0;
    if (j>_cons_j) {
        for (int l=_cons_j; l<j; ++l) {
            additionalCols += pow2i<T>(l);
        }
        additionalRows = 2*additionalCols;
    } else if (j<_cons_j) {
        for (int l=_cons_j-1; l>=j; --l) {
            additionalCols -= pow2i<T>(l);
        }
        additionalRows = 2*additionalCols;
    }

    assert(col>=_firstCol);
    assert(col<=_lastCol+additionalCols);

    if (col<=left.lastIndex()) {
        return left(col);
    }

    if (col>=right.firstIndex()+additionalCols) {
        const DenseVector<Array<T> > &rightCol = right(col-additionalCols);
        return rightCol( _ , rightCol.firstIndex()+additionalRows);
    }

    return (col>(_firstCol+_lastCol+additionalCols)/2) ?
          rightband( _ ,leftband.firstIndex() + 2*(col-left.lastIndex()-1))
        : leftband( _ , leftband.firstIndex() + 2*(col-left.lastIndex()-1));
}

template <typename T, Construction Cons>
Range<int>
RefinementMatrix<T,Interval,Cons>::rows() const
{
    return _(firstRow(), lastRow());
}

template <typename T, Construction Cons>
Range<int>
RefinementMatrix<T,Interval,Cons>::cols() const
{
    return _(firstCol(), lastCol());
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::numRows() const
{
    return lastRow()-firstRow()+1;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::numCols() const
{
    return lastCol()-firstCol()+1;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::firstRow() const
{
    return _firstRow;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::lastRow() const
{
    return _lastRow + _additionalRows;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::firstCol() const
{
    return _firstCol;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::lastCol() const
{
    return _lastCol + _additionalCols;
}

template <typename T, Construction Cons>
int
RefinementMatrix<T,Interval,Cons>::level() const
{
    return _j;
}

// TODO: consider setLevel as private(!) friend method or mra/basis.
template <typename T, Construction Cons>
void
RefinementMatrix<T,Interval,Cons>::setLevel(int j) const
{
    if (j<_j) {
        assert(j>=min_j0);
        for (int l=_j-1; l>=j; --l) {
            _additionalCols -= pow2i<T>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
    if (j>_j) {
        for (int l=_j; l<j; ++l) {
            _additionalCols += pow2i<T>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
}

template <typename T, Construction Cons>
void
RefinementMatrix<T,Interval,Cons>::_extractMasks(
                                    const GeMatrix<FullStorage<T,ColMajor> > &A)
{
    // extract left block
    for (int c=A.firstCol(); c<A.firstCol()+left.length(); ++c) {
        int r = A.lastRow();
        while (fabs(A(r,c))<=1e-12) {
            --r;
            assert(r>=A.firstRow());
        }
        left(c) = A(_(A.firstRow(),r),c);
        left(c).engine().changeIndexBase(A.firstRow());
    }
    
    // extract right block
    for (int c=A.lastCol()-right.length()+1; c<=A.lastCol(); ++c) {
        int r = A.firstRow();
        while (fabs(A(r,c))<=1e-12) {
            ++r;
            assert(r<=A.lastRow());
        }
        right(c) = A(_(r,A.lastRow()), c);
        right(c).engine().changeIndexBase(r);
    }
    
    // extract band (left to middle)
    int c = A.firstCol()+left.length();
    int first = A.firstRow();
    while (fabs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }    
    int last = A.lastRow();
    while (fabs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    leftband = A(_(first,last), c);
    leftband.engine().changeIndexBase(first);
#ifdef CHECK_INTERVAL_CONSTRUCTION
    for (++c; c<=(A.firstCol()+A.lastCol())/2; ++c) {
        int i=leftband.firstIndex();
        first += 2; last += 2;
        for (int r=first; r<=last; ++r, ++i) {
            assert(fabs(leftband(i)-A(r,c))<=1e-12);
        }
    }
#endif

    // extract band (middle to right)
    c = A.lastCol()-right.length();
    first = A.firstRow();
    while (fabs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }    
    last = A.lastRow();
    while (fabs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    rightband = A(_(first,last), c);
    rightband.engine().changeIndexBase(first);
#ifdef CHECK_INTERVAL_CONSTRUCTION
    for (--c; c>(A.firstCol()+A.lastCol())/2; --c) {
        int i=rightband.firstIndex();
        first -= 2; last -= 2;
        for (int r=first; r<=last; ++r, ++i) {
            assert(fabs(rightband(i)-A(r,c))<=1e-12);
        }
    }
    assert(leftband.length()==rightband.length());
//    assert(leftband.firstIndex()-A.firstRow()==rightband.lastIndex()-A.lastRow());
#endif

}

//------------------------------------------------------------------------------

template <typename X, Construction Cons, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,Interval,Cons> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y)
{
    typedef typename X::ElementType T;
    assert(alpha==T(1));
    assert(x.engine().stride()==1);
    assert(y.engine().stride()==1);

    if (transA==NoTrans) {
        assert(A.numCols()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.rows()) || y.engine().fill(T(0));
        } else {
            assert(y.length()==A.numRows());
            y.engine().changeIndexBase(A.firstRow());
        }

        // left upper block
        int ix = x.firstIndex();
        for (int c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++ix) {
            int n = A.left(c).length();
            cxxblas::axpy(n, 
                          x(ix),
                          A.left(c).engine().data(), 1,
                          y.engine().data(), 1);
        }

        // central band (up to middle)
        int iy = A.leftband.firstIndex()-A.firstRow();
        int n = A.leftband.length();
        int middle = iceil<int>(x.length()/2.);
        for (int c=A.left.lastIndex()+1; c<=middle; ++c, iy+=2, ++ix) {
            cxxblas::axpy(n,
                          x(ix),
                          A.leftband.engine().data(), 1,
                          y.engine().data()+iy, 1);
        }
        // central band (right of middle)
        int end = A.left.firstIndex() + x.length() - A.right.length();
        for (int c=middle+1; c<end; ++c, iy+=2, ++ix) {
            cxxblas::axpy(n,
                          x(ix),
                          A.rightband.engine().data(), 1,
                          y.engine().data()+iy, 1);
        }

        // right lower block
        for (int c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++ix) {
            int n = A.right(c).length();
            cxxblas::axpy(n, 
                          x(ix),
                          A.right(c).engine().data(), 1,
                          y.engine().data()+y.length()-1-n+1, 1);
        }
    } else { // transA==Trans
        assert(A.numRows()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.cols());
        } else {
            assert(y.length()==A.numCols());
            y.engine().changeIndexBase(A.firstCol());
        }
        int iy = y.firstIndex();
        // left upper block
        for (int c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++iy) {
            int n = A.left(c).length();
            cxxblas::dot(n,
                         A.left(c).engine().data(), 1,
                         x.engine().data(), 1, 
                         y(iy));
        }

        // central band (up to middle)
        int middle = y.length()/2;
        int ix = A.leftband.firstIndex() - A.firstRow();
        for (int i=A.left.length()+1; i<=middle; ++i, ix+=2, ++iy) {
            cxxblas::dot(A.leftband.length(),
                         A.leftband.engine().data(), 1,
                         x.engine().data()+ix, 1,
                         y(iy));
        }
        // central band (right of middle)
        int end = y.length() - A.right.length();
        for (int i=middle+1; i<=end; ++i, ix+=2, ++iy) {
            cxxblas::dot(A.rightband.length(),
                         A.rightband.engine().data(), 1,
                         x.engine().data()+ix, 1,
                         y(iy));
        }
        // right lower block
        for (int c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++iy) {
            int n = A.right(c).length();
            cxxblas::dot(n, 
                         A.right(c).engine().data(), 1,
                         x.engine().data() + A.numRows() - n, 1, 
                         y(iy));
        }
    }
}

} // namespace flens

