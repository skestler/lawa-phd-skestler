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

namespace cxxblas {

#ifdef HAVE_CBLAS
// srot
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
rot(IndexType n, float *x, IndexType incX, float *y, IndexType incY,
    float c, float s)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_srot");

    cblas_srot(n, x, incX, y, incY, c, s);
}

// drot
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
rot(IndexType n, double *x, IndexType incX, double *y, IndexType incY,
    double c, double s)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_drot");

    cblas_drot(n, x, incX, y, incY, c, s);
}

// srotg
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
rotg(float &a, float &b, float &c, float &s)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_srotg");

    cblas_srotg(&a, &b, &c, &s);
}

// drotg
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
rotg(double &a, double &b, double &c, double &s)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_drotg");

    cblas_drotg(&a, &b, &c, &s);
}

#endif // HAVE_CBLAS

} // namespace cxxblas
