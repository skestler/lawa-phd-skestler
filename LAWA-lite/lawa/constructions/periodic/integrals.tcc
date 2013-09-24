/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#include <lawa/math/math.h>
#include <lawa/constructions/support.h>

namespace lawa {

// reduce periodic integral to equivalent integral on R.
template <typename T>
void
_adapt_k(const Support<T> &supp1,
         const Support<T> &supp2, 
         int j1, int &k1, int j2, int &k2)
{
    if (supp1.l2>1.) {
        if (supp2.l2<supp1.l1) {
        k1 -= pow2i<T>(j1);
        }
    } else if (supp2.l2>1.) {
        if (supp1.l2<supp2.l1) {
            k2 -= pow2i<T>(j2);
        }
    } else if (supp1.l1<0.) {
        if (supp2.l1>supp1.l2) {
            k1 += pow2i<T>(j1);
        }
    } else if (supp2.l1<0) {
        if (supp1.l1>supp2.l2) {
            k2 += pow2i<T>(j2);
        }
    }
}

} // namespace lawa

