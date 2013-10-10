/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BASIS_H
#define LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Primal,R,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = R;
        static const Construction Cons = SparseMulti;

        typedef BasisFunction<T,Primal,R,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,R,SparseMulti>       BSplineType;
        typedef Wavelet<T,Primal,R,SparseMulti>       WaveletType;

        Basis(const int d, const int j=-1);

        int
        level() const;

        void
        setLevel(int j) const;

        const BasisFunctionType &
        generator(XType xtype) const;

        MRA<T,Primal,R,SparseMulti> mra;

        const int d, j0;

        Wavelet<T,Primal,R,SparseMulti> psi;


    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/realline/sparsemulti/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BASIS_H

