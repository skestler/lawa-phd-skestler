/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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
 
namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
TensorBasis2D<SparseGrid, FirstBasis, SecondBasis>::TensorBasis2D(const FirstBasis &_basis1, 
                                                                  const SecondBasis &_basis2)  
    : first(_basis1), second(_basis2)
{
}

template<typename FirstBasis, typename SecondBasis>
int
TensorBasis2D<SparseGrid, FirstBasis, SecondBasis>::dim(const int J_x, const int J_y) const
{
    int d = first.mra.cardI(first.j0) * second.mra.cardI(J2_max(J_x, J_y, first.j0 - 1));
    for(int jx = first.j0; jx <= J1_max(J_x, J_y, second.j0-1) - 1; ++jx){
        d += first.cardJ(jx)*second.mra.cardI(second.j0);
        for(int jy = second.j0; jy <= J2_max(J_x,J_y, jx) - 1; ++jy){
            d += first.cardJ(jx) * second.cardJ(jy);
        }
    }
    return d;
}

/* Calculate maximal level in dimension 1 if level in dimension 2 is jy.
 * 
 *  Convention: Level of scaling function space (S_j) is set to j-1 
 *
 *  Isotropic Spaces: Full Diagonals
 *      C = C_x = (J_x - 1) + (j0_y - 1) = C_y (= (J_y - 1) + (j0_x - 1))
 *  Anisotropic Spaces: always include spaces with maximal levels J_x, J_y, 
 *      then fill diagonals as much as possible
 *      C = max(C_x, C_y)
 * 
 *  Constraints: j_x < J_x & j_x <= C - j_y
 *  => j_x < min(J_x, C - j_y + 1)
 */
template<typename FirstBasis, typename SecondBasis>
int 
TensorBasis2D<SparseGrid, FirstBasis, SecondBasis>::J1_max(const int J_x, const int J_y, const int jy) const
{
    return std::min(J_x, std::max(J_x + second.j0 - 2, J_y + first.j0 - 2) - jy + 1);
} 

template<typename FirstBasis, typename SecondBasis>
int 
TensorBasis2D<SparseGrid, FirstBasis, SecondBasis>::J2_max(const int J_x, const int J_y, const int jx) const
{
    return std::min(J_y, std::max(J_x + second.j0 - 2, J_y + first.j0 - 2)- jx  + 1);
}
    
} // namespace lawa
