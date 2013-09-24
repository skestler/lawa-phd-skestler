/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

namespace lawa {
    
template<typename Basis>
UniformIndex2D<Basis>::UniformIndex2D(const Basis& _basis, const int _J_x, const int _J_y)
    : basis(_basis), J_x(_J_x), J_y(_J_y),
      offsetIx(_basis.first.mra.rangeI(_basis.first.j0).firstIndex() - 1), 
      offsetIy(_basis.second.mra.rangeI(_basis.second.j0).firstIndex() - 1), 
      offsetJx(_basis.first.rangeJ(_basis.first.j0).firstIndex() - 1),
      offsetJy(_basis.second.rangeJ(_basis.second.j0).firstIndex() - 1)
{
}

template<typename Basis>
int
UniformIndex2D<Basis>::operator()(XType xtype_x, int jx, int kx,
                                  XType xtype_y, int jy, int ky) const
{

    if(xtype_x == XBSpline){
         if(xtype_y == XBSpline){
             return (kx-offsetIx-1)*basis.second.mra.cardI(J_y)
                   + ky-offsetIy;
         }
         else{
             return  (kx-offsetIx-1)*basis.second.mra.cardI(J_y)
                     +  basis.second.mra.cardI(jy)
                     +  ky-offsetJy;
         }
     }
     else{
         if(xtype_y == XBSpline){
             return  basis.dim(jx, J_y)
                     + (kx-offsetJx-1)*basis.second.mra.cardI(J_y)
                     + ky-offsetIy;

         }
         else{
             return  basis.dim(jx, J_y)
                     + (kx-offsetJx-1)*basis.second.mra.cardI(J_y)
                     +  basis.second.mra.cardI(jy)
                     +  ky-offsetJy;
         }
     }

    /*
    if(xtype_x == XBSpline){
        if(xtype_y == XBSpline){
            return (kx-offsetIx-1)*basis.second.mra.cardI(basis.second.j0)
                  + ky-offsetIy;
        }
        else{
            return  basis.first.mra.cardI(basis.first.j0) * basis.second.mra.cardI(jy)
                    + (kx-offsetIx-1)*basis.second.cardJ(jy)
                    +  ky-offsetJy;
        }
    }
    else{
        if(xtype_y == XBSpline){
            return  basis.dim(jx, J_y)
                    + (kx-offsetJx-1)*basis.second.mra.cardI(basis.second.j0)
                    +  ky-offsetIy;

        }
        else{
            return  basis.dim(jx, J_y)
                    + basis.first.cardJ(jx) * basis.second.mra.cardI(jy)
                    + (kx-offsetJx-1)*basis.second.cardJ(jy)
                    +  ky-offsetJy;
        }
    }
    */

}
    
} // namespace lawa

