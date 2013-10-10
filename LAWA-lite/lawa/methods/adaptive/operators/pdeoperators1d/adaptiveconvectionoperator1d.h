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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVECONVECTIONOPERATOR1D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVECONVECTIONOPERATOR1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>
#include <lawa/preconditioners/preconditioners.h>
//#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct AdaptiveConvectionOperator1D
{
    typedef Basis<T,Side,Domain,Cons>                                Basis1D;
    typedef CompressionPDE1D<T, Basis1D>                             Compression1D;
    typedef NoPreconditioner<T,Index1D>                              NoPreconditioner1D;
    typedef ConvectionOperator1D<T, Basis1D>                         ConvectionOp1D;
    typedef MapMatrix<T, Index1D, ConvectionOp1D,
                     Compression1D, NoPreconditioner1D>              DataConvectionOp1D;

    AdaptiveConvectionOperator1D(const Basis1D &_basis1d);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    void
    clear();

    const Basis1D               &basis1d;
    Compression1D               compression1d;
    const ConvectionOp1D        convection_op1d;
    NoPreconditioner1D          prec1d;
    DataConvectionOp1D          data;

};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveconvectionoperator1d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVECONVECTIONOPERATOR1D_H

