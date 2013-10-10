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

#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_OPTIMIZED_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_OPTIMIZED_H 1

#include <iostream>
#include <vector>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/algorithms/linearsystemsolvers.h>
#include <lawa/methods/adaptive/postprocessing/postprocessing.h>

namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
class S_ADWAV_Optimized {

    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;

    public:
        S_ADWAV_Optimized(AdaptiveOperator &A, RHS &F, PP_AdaptiveOperator &PP_A, PP_RHS &PP_F,
                          T contraction, T start_threshTol,
                          T _linTol, T _resTol, int _NumOfIterations, int _MaxItsPerThreshTol,
                          T _eps=1e-2);

        //solver for stationary problems
        void
        solve(const IndexSet<Index> &Initial_Lambda, Coefficients<Lexicographical,T,Index> &u,
              const char *linsolvertype, const char *filename, int assemble_matrix, T H1norm=0.);

        std::vector<T>               residuals;
        std::vector<T>               times;
        std::vector<T>               linsolve_iterations;
        std::vector<T>               toliters;

    private:
        AdaptiveOperator &A;
        RHS &F;
        PP_AdaptiveOperator &PP_A;
        PP_RHS &PP_F;
        T contraction, threshTol, linTol, resTol;
        int NumOfIterations;
        int MaxItsPerThreshTol;
        T eps;

};

} //namespace lawa

#include <lawa/methods/adaptive/solvers/s_adwav_optimized.tcc>

#endif //LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_OPTIMIZED_H

