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

namespace flens {

template <typename T>
GeMatrix<FullStorage<T,cxxblas::ColMajor> >
inv(const GeMatrix<FullStorage<T,cxxblas::ColMajor> > &A)
{
    assert(A.numRows()==A.numCols());
    
    // inversion using QR ...
    GeMatrix<FullStorage<T,cxxblas::ColMajor> > I(A.numRows(),A.numCols()), 
                                                Tmp1 = A, Tmp2, InvA;
    I.diag(0) = 1;
    flens::DenseVector<Array<T> > tau;
    qrf(Tmp1, tau);
    Tmp2 = Tmp1;
    orgqr(Tmp1, tau);

    blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,Tmp1,I,0.,InvA);
    blas::sm(cxxblas::Left,cxxblas::NoTrans,1.,Tmp2.upper(),InvA);

    return InvA;
}

} // namespace flens

