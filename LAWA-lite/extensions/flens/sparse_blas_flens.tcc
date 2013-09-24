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

#include <cassert>
#include <extensions/flens/sparse_blas.h>

namespace flens { namespace blas {

template <typename T, typename VX, typename VY>
void
mv(cxxblas::Transpose trans, T alpha, const SparseGeMatrix<CRS<T> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::ElementType beta, DenseVector<VY> &y)
{
    assert(ADDRESS(y)!=ADDRESS(x));
    if (trans==cxxblas::NoTrans) {
        assert(A.numCols()==x.length());
    }
    if (trans==cxxblas::Trans) {
        assert(A.numRows()==x.length());
    }
    int yLength = (trans==cxxblas::NoTrans) ? A.numRows()
                                            : A.numCols();

    assert((beta==0) || (y.length()==yLength));

    if (y.length()!=yLength) {
        y.engine().resize(yLength);
    }

    assert(x.engine().stride()==1);
    assert(y.engine().stride()==1);

    crs_gemv(trans, A.numRows(), A.numCols(),
             alpha,
             A.engine().values.engine().data(),
             A.engine().rows.engine().data(),
             A.engine().columns.engine().data(),
             x.engine().data(),
             beta, y.engine().data());
}

// sparse_symv
template <typename T, CRS_Storage Storage, typename VX, typename VY>
void
mv(T alpha, const SparseSyMatrix<CRS<T, Storage> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::ElementType beta, DenseVector<VY> &y)
{
    assert(A.dim()==x.length());
    assert((beta==0) || (A.dim()==y.length()));

    if (y.length()!=A.dim()) {
        y.engine().resize(A.dim());
    }

    assert(x.engine().stride()==1);
    assert(y.engine().stride()==1);

    if (Storage==CRS_General) {
        crs_gemv(cxxblas::NoTrans, A.dim(), A.dim(),
                 alpha,
                 A.engine().values.data(),
                 A.engine().rows.data(),
                 A.engine().columns.data(),
                 x.data(),
                 beta, y.data());
    } else {
        cxxblas::StorageUpLo upLo = (Storage==CRS_UpperTriangular) ? cxxblas::Upper : cxxblas::Lower;
        crs_symv(upLo, A.dim(), alpha,
                 A.engine().values.data(),
                 A.engine().rows.data(),
                 A.engine().columns.data(),
                 x.data(),
                 beta, y.data());
    }
}

template <typename T>
void
my_mm_t(const SparseGeMatrix<CRS<T> > &A,
        SparseGeMatrix<CRS<T> > &B)
{
    assert(B.numRows() == A.numCols());
    assert(B.numCols() == A.numCols());
    int m = A.numRows();
    const T *a;
    const int *ia, *ja;
    a  = A.engine().values.engine().data();
    ia = A.engine().rows.engine().data();
    ja = A.engine().columns.engine().data();
    a = a-1;
    ia = ia-1;
    ja = ja-1;

    for (int i_t=1; i_t<=m; ++i_t) {
        for (int k_t=ia[i_t]; k_t<ia[i_t+1]; ++k_t) {
            //i_t = row index of A = column index of A^t
            //k_t = column index of A = row index of A^t
            T a_k_t_i_t = a[k_t];
            for (int i=1; i<=m; ++i) {
                for (int k=ia[i]; k<ia[i+1]; ++k) {
                    //i = row index of A
                    //k = column index of A
                    if (i_t == i) B(ja[k_t],ja[k]) = a_k_t_i_t * a[k];
                }
            }
        }
    }
    B.finalize();
}

template <typename T>
void
my_mm_At_A(const SparseGeMatrix<CRS<T> > &A,
           SparseGeMatrix<CRS<T> > &B)
{
    assert(B.numRows() == A.numRows());
    assert(B.numCols() == A.numRows());
    int m = A.numRows();
    const T *a;
    const int *ia, *ja;
    a  = A.engine().values.engine().data();
    ia = A.engine().rows.engine().data();
    ja = A.engine().columns.engine().data();
    a = a-1;
    ia = ia-1;
    ja = ja-1;

    for (int i_t=1; i_t<=m; ++i_t) {
        //int col_At_left =  ja[ia[i_t]];
        //int col_At_right = ja[ia[i_t+1]-1];
        for (int i=1; i<=i_t; ++i) {
            T a_i_t_i = 0.;
            //int row_A_left  =  ja[ia[i]];
            //int row_A_right  = ja[ia[i+1]-1];
            //int left  = std::max(col_At_left,row_A_left);
            //int right = std::min(col_At_right,row_A_right);
            for (int k_t=ia[i_t]; k_t<=ia[i_t+1]-1; ++k_t) {
                //if (ja[k_t]<left)  continue;
                //if (ja[k_t]>right) break;
                T a_k_t = a[k_t];
                for (int k=ia[i]; k<=ia[i+1]-1; ++k) {
                    //if (ja[k]<left)  continue;
                    //if (ja[k]>right) break;
                    if (ja[k_t] == ja[k]) a_i_t_i += a_k_t * a[k];
                }
            }
            B(i_t,i) = a_i_t_i;
            if (i != i_t) B(i,i_t) = a_i_t_i;
        }
    }
    B.finalize();
}

// sparse_gemm (matrix B dense!)
template <typename T, CRS_Storage Storage>
void
mm(cxxblas::Transpose transA,
   cxxblas::Transpose transB,
   T alpha,
   const SparseGeMatrix<CRS<T, Storage> > &A,
   const GeMatrix<FullStorage<T, cxxblas::ColMajor> > &B,
   T beta, GeMatrix<FullStorage<T, cxxblas::ColMajor> > &C)
{
    assert(Storage==CRS_General);
    assert(transB==cxxblas::NoTrans); // Transposition not provided for B.
    int m = (transA==cxxblas::NoTrans) ? A.numRows() : A.numCols();
    int n = B.numCols();

    assert((beta==0) || (C.numRows()==m));
    assert((beta==0) || (C.numCols()==n));

    C.engine().resize(m,n);

    int r = A.numRows();
    DenseVector<Array<int> > pointerE(r,1);
    pointerE(Range<int>(1,r-1)) = A.engine().rows(Range<int>(2,r));
    pointerE(r) = A.engine().values.length()+1;

    char matdescra[] = {'G', 'U', 'N', 'F', ' ', ' '};

    csrmm(transA, A.numRows(), C.numCols(), A.numCols(), alpha, matdescra,
          A.engine().values.engine().data(), A.engine().columns.engine().data(),
          A.engine().rows.engine().data(), pointerE.engine().data(),
          B.engine().data(), B.leadingDimension(), 
          beta, C.engine().data(), C.leadingDimension());
}

} } // namespase blas, flens

