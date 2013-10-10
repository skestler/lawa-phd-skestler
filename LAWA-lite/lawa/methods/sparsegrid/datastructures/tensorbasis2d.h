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

#ifndef LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORBASIS2D_H
#define LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORBASIS2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
struct TensorBasis2D<SparseGrid, FirstBasis, SecondBasis>
{
    typedef typename FirstBasis::T T;
    typedef FirstBasis FirstBasisType;
    typedef SecondBasis SecondBasisType;
    
    TensorBasis2D(const FirstBasis &_basis1, const SecondBasis &_basis2);

    const FirstBasis &first;
    const SecondBasis &second;
    
    int
    dim(const int J_x, const  int J_y) const;
    
    /* J1_max and J2_max return the maximal level in dimension 1 (2) 
     * if the level in the other dimension is jy (jx).
     * This allows to use the uniform assembling code also for sparsegrid
     * tensorbases.  
     * Here: only "diagonal" spaces (i.e. jx + jy <= C)
     *  For anisotropic bases, we fill the diagonals "as much as possible", i.e.
     *  C = max(C_x, C_y).
     */
    int 
    J1_max(const int J_x, const int J_y, const int jy) const;
    
    int 
    J2_max(const int J_x, const int J_y, const int jx) const;
};
    
    
} // namespace lawa

#include <lawa/methods/sparsegrid/datastructures/tensorbasis2d.tcc>


#endif // LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORBASIS2D_H
