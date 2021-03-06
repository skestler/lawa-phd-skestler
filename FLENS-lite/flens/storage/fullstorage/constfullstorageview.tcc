/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

namespace flens {

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(
                                        const ElementType *data,
                                        const Allocator &allocator,
                                        IndexType numRows, IndexType numCols,
                                        IndexType leadingDimension,
                                        IndexType firstRow,
                                        IndexType firstCol)
    : _data(data), _allocator(allocator), _numRows(numRows), _numCols(numCols),
      _leadingDimension(leadingDimension), _firstRow(0), _firstCol(0)
{
    changeIndexBase(firstRow, firstCol);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(
                                                const ConstFullStorageView &rhs)
    : _data(rhs._data),
      _allocator(rhs._allocator),
      _numRows(rhs._numRows), _numCols(rhs._numCols),
      _leadingDimension(rhs._leadingDimension),
      _firstRow(rhs._firstRow), _firstCol(rhs._firstCol)
{
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
template <typename RHS>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(const RHS &rhs)
    : _data(rhs.data()),
      _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _leadingDimension(rhs.leadingDimension()),
      _firstRow(0), _firstCol(0)
{
    changeIndexBase(rhs.firstRow(), rhs.firstCol());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::~ConstFullStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ElementType &
ConstFullStorageView<T, Order, I, A>::operator()(IndexType row,
                                                 IndexType col) const
{
    ASSERT(row>=_firstRow);
    ASSERT(row<_firstRow+_numRows);
    ASSERT(col>=_firstCol);
    ASSERT(col<_firstCol+_numCols);

    if (Order==cxxblas::ColMajor) {
        return _data[col*_leadingDimension+row];
    }
    return _data[row*_leadingDimension+col];
}

//-- methods -------------------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::firstRow() const
{
    return _firstRow;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::firstCol() const
{
    return _firstCol;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::numRows() const
{
    return _numRows;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::numCols() const
{
    return _numCols;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::leadingDimension() const
{
    return _leadingDimension;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::strideRow() const
{
    return (Order==cxxblas::ColMajor) ? 1
                                      : leadingDimension();
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::strideCol() const
{
    return (Order==cxxblas::ColMajor) ? leadingDimension()
                                      : 1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ElementType *
ConstFullStorageView<T, Order, I, A>::data() const
{
    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::Allocator &
ConstFullStorageView<T, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
ConstFullStorageView<T, Order, I, A>::changeIndexBase(IndexType firstRow,
                                                      IndexType firstCol)
{
    if (Order==cxxblas::RowMajor) {
        _data = data() - (firstRow*leadingDimension() + firstCol);
    }
    if (Order==cxxblas::ColMajor) {
        _data = data() - (firstCol*leadingDimension() + firstRow);
    }
    _firstRow = firstRow;
    _firstCol = firstCol;
}

// view of rectangular part
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const ConstFullStorageView<T, Order, I, A>
ConstFullStorageView<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                           IndexType toRow, IndexType toCol,
                                           IndexType firstViewRow,
                                           IndexType firstViewCol) const
{
    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    return ConstView(&(this->operator()(fromRow, fromCol)),// data
                     allocator(),                          // allocator
                     toRow-fromRow+1,                      // # rows
                     toCol-fromCol+1,                      // # cols
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol);                        // firstCol
}

// view of single row
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                              IndexType firstViewIndex) const
{
    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(&(this->operator()(row, _firstCol))-firstViewIndex,
                          allocator(),
                          numCols(),
                          strideCol(),
                          firstViewIndex);
}

// view of single columns
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewCol(IndexType col,
                                              IndexType firstViewIndex) const
{
    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(&(this->operator()(_firstRow, col))-firstViewIndex,
                          allocator(),
                          numRows(),
                          strideRow(),
                          firstViewIndex);
}

// view of d-th diagonal
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewDiag(IndexType d,
                                               IndexType firstViewIndex) const
{
    IndexType col = firstCol() + ( (d>0) ? d : 0 );
    IndexType row = firstRow() + ( (d>0) ? 0 : -d );
    return ConstArrayView(&(this->operator()(row,col)),
                          allocator(),
                          std::min(numRows(),numCols()) - std::abs(d),
                          leadingDimension()+1,
                          firstViewIndex);
}

} // namespace flens
