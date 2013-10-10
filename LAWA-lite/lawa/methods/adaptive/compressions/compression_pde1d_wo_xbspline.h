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


#ifndef  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE1D_WO_XBSPLINE_H
#define  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE1D_WO_XBSPLINE_H 1

#include <lawa/settings/settings.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/aux/timer.h>

namespace lawa {

//only for R-Basis!!
template <typename T, typename Basis>
struct CompressionPDE1D_WO_XBSpline
{
    const Basis &basis;
    short s_tilde, jmin, jmax;

    CompressionPDE1D_WO_XBSpline(const Basis &_basis);

    void
    setParameters(const IndexSet<Index1D> &LambdaRow);

    IndexSet<Index1D>
    SparsityPattern(const Index1D &lambda, const IndexSet<Index1D> &Lambda, int J=-1);

};

} // namespace lawa

#include <lawa/methods/adaptive/compressions/compression_pde1d_wo_xbspline.tcc>

#endif //  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE1D_WO_XBSPLINE_H

