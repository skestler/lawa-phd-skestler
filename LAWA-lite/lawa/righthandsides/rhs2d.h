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

#ifndef LAWA_RIGHTHANDSIDES_RHS2D_H
#define LAWA_RIGHTHANDSIDES_RHS2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {
    
template <typename T>
struct Rhs2D {
	
    virtual T
    operator()(const Index2D &index) const = 0;
    
    virtual T
    operator()(XType xtype_x, int j_x, long k_x,
               XType xtype_y, int j_y, long k_y) const = 0;
};
    
} // namespace lawa

#endif // LAWA_RIGHTHANDSIDES_RHS2D_H
