/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Mario Rometsch.

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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_TCC 1

namespace lawa {

template <typename T>
T
_constant_bspline_inner_evaluator0(T x, unsigned short deriv) {

    if(deriv == 0){
        if(0. <= x && x <1.){
            return 1.L;
        } else {
            return 0.L;
        }
    }
    else {
        return 0.L;
    }
}

// Attention: do not evaluate at x=1! Refinement of wavelet by scaling is not possible at this point
// since if we evaluate the scaling function at x=1, then function values double at interface when
// calculating phi(2x) + phi(2x-1)!
template <typename T>
T
_constant_wavelet_inner_evaluator0(T x, unsigned short deriv) {

    if(deriv == 0){
        if(0. <= x && x <0.5){
            return 1.L;
        } else if(0.5 <= x && x <1.) {
            return -1.L;
        }
        else return 0.L;
    }
    else {
        return 0.L;
    }
}

}   // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI__CONSTANT_EVALUATOR_TCC
