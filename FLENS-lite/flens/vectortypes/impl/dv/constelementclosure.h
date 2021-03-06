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

#ifndef FLENS_VECTORTYPES_IMPL_DV_CONSTELEMENTCLOSURE_H
#define FLENS_VECTORTYPES_IMPL_DV_CONSTELEMENTCLOSURE_H 1

#include <flens/aux/constref.h>
#include <flens/scalartypes/scalar.h>

namespace flens { namespace densevector {

template <typename V, typename I = typename V::IndexVariable>
class ConstElementClosure
    : public Scalar<ConstElementClosure<V, I> >
{
    public:
        typedef V                               Vector;
        typedef typename Vector::ElementType    ElementType;
        typedef I                               IndexVariable;

        ConstElementClosure(const Vector &vector,
                            const IndexVariable &index);

        const ElementType &
        value() const;

    private:
        const Vector         &_vector;
        const IndexVariable  &_index;
        int _id;
};

} } // namespace densevector, flens

#include <flens/vectortypes/impl/dv/constelementclosure.tcc>

#endif // FLENS_VECTORTYPES_IMPL_DV_CONSTELEMENTCLOSURE_H
