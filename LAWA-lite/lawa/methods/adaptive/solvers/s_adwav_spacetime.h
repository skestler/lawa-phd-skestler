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
#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_SPACETIME_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_SPACETIME_H 1

#include <iostream>
#include <vector>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/postprocessing/postprocessing.h>

namespace lawa {

template <typename T, typename Index, typename SpaceIndex, typename Basis, typename MA, typename RHSOperator, typename RHSInitialCond>
class S_ADWAV_SPACETIME {
    public:
        S_ADWAV_SPACETIME(const Basis &basis, MA &A, RHSOperator &F, RHSInitialCond &U0, T contraction, T _threshTol,
                          T _linTol, T _resTol=1e-4, int _NumOfIterations=10, int MaxItsPerThreshTol=5, T eps=1e-2);

        void solve_cgls(const IndexSet<Index> &Initial_Lambda);

        std::vector<Coefficients<Lexicographical,T,Index> > solutions;
        std::vector<T>               residuals;
        std::vector<T>               times;

    private:
        const Basis &basis;
        MA &A;
        RHSOperator &F;
        RHSInitialCond &U0;
        T contraction, threshTol, linTol, resTol;
        int NumOfIterations;
        int MaxItsPerThreshTol;
        T eps;

};

} //namespace lawa

#include <lawa/methods/adaptive/solvers/s_adwav_spacetime.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_SPACETIME_H

