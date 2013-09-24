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

#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

void
getMinAndMaxLevel(const IndexSet<Index1D> &Lambda, int &jmin, int &jmax)
{
    typedef IndexSet<Index1D>::const_iterator set1d_const_it;
    set1d_const_it it = Lambda.begin();
    jmin = (*it).j;
    jmax = (*it).j;
    for (set1d_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        jmin = std::min(int((*lambda).j),jmin);
        jmax = std::max(int((*lambda).j),jmax);
    }
}

void
split(const IndexSet<Index2D> &Lambda, IndexSet<Index1D> &Lambda_x, IndexSet<Index1D> &Lambda_y)
{
    typedef IndexSet<Index2D>::const_iterator set2d_const_it;
    for (set2d_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        Lambda_x.insert((*lambda).index1);
        Lambda_y.insert((*lambda).index2);
    }
}

IndexSet<Index1D>
extractSpaceIndices(const IndexSet<Index2D> &Lambda)
{
    typedef IndexSet<Index2D>::const_iterator set2d_const_it;
    IndexSet<Index1D> ret;
    for (set2d_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        ret.insert((*lambda).index2);
    }
    return ret;
}
    
} // namespace lawa
