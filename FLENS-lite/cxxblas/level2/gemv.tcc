/*
 *   Copyright (c) 2009, Michael Lehn
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

#include <complex>
#include <cxxblas/level1/level1.h>

namespace cxxblas {

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
gemv_generic(StorageOrder order, Transpose trans,
             IndexType m, IndexType n,
             const ALPHA &alpha,
             const MA *A, IndexType ldA,
             const VX *x, IndexType incX,
             const BETA &beta,
             VY *y, IndexType incY)
{
    if (order==ColMajor) {
        trans = Transpose(trans^Trans);
        gemv_generic(RowMajor, trans, n, m, alpha, A, ldA,
                     x, incX, beta, y, incY);
        return;
    }
    if ((trans==NoTrans) || (trans==Conj)) {
        scal(m, beta, y, incY);
        if (trans==Conj) {
            for (IndexType i=0, iY=0; i<m; ++i, iY+=incY) {
                VY _y;
                dot(n, A+i*ldA, IndexType(1), x, incX, _y);
                y[iY] += alpha*_y;
            }
        } else {
            for (IndexType i=0, iY=0; i<m; ++i, iY+=incY) {
                VY _y;
                dotu(n, A+i*ldA, IndexType(1), x, incX, _y);
                y[iY] += alpha*_y;
            }
        }
    } else {
        scal(n, beta, y, incY);
        if (trans==ConjTrans) {
            for (IndexType i=0, iY=0; i<n; ++i, iY+=incY) {
                VY _y;
                dot(n, A+i, ldA, x, incX, _y);
                y[iY] += alpha*_y;
            }
        } else {
            for (IndexType i=0, iY=0; i<n; ++i, iY+=incY) {
                VY _y;
                dotu(n, A+i, ldA, x, incX, _y);
                y[iY] += alpha*_y;
            }
        }
    }
}

//------------------------------------------------------------------------------

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
gemv(StorageOrder order, Transpose trans,
     IndexType m, IndexType n,
     const ALPHA &alpha,
     const MA *A, IndexType ldA,
     const VX *x, IndexType incX,
     const BETA &beta,
     VY *y, IndexType incY)
{
    gemv_generic(order, trans, m, n, alpha, A, ldA, x, incX, beta, y, incY);
}

#ifdef HAVE_CBLAS

// sgemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
gemv(StorageOrder order, Transpose trans,
     IndexType m, IndexType n,
     float alpha,
     const float *A, IndexType ldA,
     const float *x, IndexType incX,
     float beta,
     float *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_sgemv");

    cblas_sgemv(CBLAS::value(order), CBLAS::value(trans),
                m,  n,
                alpha,
                A, ldA,
                x, incX,
                beta,
                y, incY);
}

// dgemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
gemv(StorageOrder order, Transpose trans,
     IndexType m, IndexType n,
     double alpha,
     const double *A, IndexType ldA,
     const double *x, IndexType incX,
     double beta,
     double *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_dgemv");

    cblas_dgemv(CBLAS::value(order), CBLAS::value(trans),
                m,  n,
                alpha,
                A, ldA,
                x, incX,
                beta,
                y, incY);
}

// cgemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
gemv(StorageOrder order, Transpose trans,
     IndexType m, IndexType n,
     ComplexFloat &alpha,
     const ComplexFloat *A, IndexType ldA,
     const ComplexFloat *x, IndexType incX,
     ComplexFloat &beta,
     ComplexFloat *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_cgemv");

    cblas_cgemv(CBLAS::value(order), CBLAS::value(trans),
                m,  n,
                reinterpret_cast<const float *>(&alpha),
                reinterpret_cast<const float *>(A), ldA,
                reinterpret_cast<const float *>(x), incX,
                reinterpret_cast<const float *>(&beta),
                reinterpret_cast<const float *>(y), incY);
}

// zgemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
gemv(StorageOrder order, Transpose trans,
     IndexType m, IndexType n,
     ComplexDouble &alpha,
     const ComplexDouble *A, IndexType ldA,
     const ComplexDouble *x, IndexType incX,
     ComplexDouble &beta,
     ComplexDouble *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_zgemv");

    cblas_zgemv(CBLAS::value(order), CBLAS::value(trans),
                m,  n,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(A), ldA,
                reinterpret_cast<const double *>(x), incX,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<const double *>(y), incY);
}

#endif // HAVE_CBLAS

} // namespace cxxblas
