/*
 *   Copyright (c) 2010, Michael Lehn
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

namespace cxxblas {

#ifdef HAVE_CBLAS

// chemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
hemv(StorageOrder order, StorageUpLo upLo,
     IndexType n, ComplexDouble &alpha,
     const ComplexFloat *A, IndexType ldA,
     const ComplexFloat *x, IndexType incX,
     ComplexDouble &beta,
     ComplexFloat *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_chemv");

    cblas_chemv(CBLAS::value(order), CBLAS::value(upLo), n, 
                reinterpret_cast<const float *>(&alpha),
                reinterpret_cast<const float *>(A), ldA,
                reinterpret_cast<const float *>(x), incX,
                reinterpret_cast<const float *>(&beta),
                reinterpret_cast<const float *>(y), incY);
}

// zhemv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
hemv(StorageOrder order, StorageUpLo upLo,
     IndexType n, ComplexDouble &alpha,
     const ComplexDouble *A, IndexType ldA,
     const ComplexDouble *x, IndexType incX,
     ComplexDouble &beta,
     ComplexDouble *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_zhemv");

    cblas_zhemv(CBLAS::value(order), CBLAS::value(upLo), n,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(A), ldA,
                reinterpret_cast<const double *>(x), incX,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<const double *>(y), incY);
}

#endif // HAVE_CBLAS

} // namespace cxxblas
