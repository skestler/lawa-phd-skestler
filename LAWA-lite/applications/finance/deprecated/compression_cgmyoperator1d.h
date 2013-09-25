/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef  APPLICATIONS_FINANCE_COMPRESSION_CGMYOPERATOR1D_H
#define  APPLICATIONS_FINANCE_COMPRESSION_CGMYOPERATOR1D_H

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/aux/timer.h>

namespace lawa {

//only for R-Basis!!
template <typename T, typename Basis>
struct CompressionCGMYOperator1D
{
};

template <typename T>
struct CompressionCGMYOperator1D<T, Basis<T,Primal,R,CDF> >
{
    typedef typename Basis<T,Primal,R,CDF>::WaveletType PrimalWavelet;

    const Basis<T,Primal,R,CDF> &basis;
    short jmin, jmax;
    T compr_c;
    T compr_alpha;
    PrimalWavelet psi;

    CompressionCGMYOperator1D(const Basis<T,Primal,R,CDF> &_basis, T Y);

    void
    setParameters(const IndexSet<Index1D> &LambdaRow);

    IndexSet<Index1D>
    SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow);

};

} // namespace lawa

#include <applications/finance/compression_cgmyoperator1d.tcc>

#endif //  APPLICATIONS_FINANCE_COMPRESSION_CGMYOPERATOR1D_H
