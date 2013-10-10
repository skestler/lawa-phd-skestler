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
#ifndef LAWA_METHODS_ADAPTIVE_POSTPROCESSING_ERROR_ESTIMATORS_H
#define LAWA_METHODS_ADAPTIVE_POSTPROCESSING_ERROR_ESTIMATORS_H 1

#include <iostream>
#include <fstream>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <vector>

namespace lawa {

//todo: Use APPLY if possible!
// Calculates an estimate for ||Au_{\Lambda}-f||_{ell2}
template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_Au_M_f(MA &A, RHS &F, const Coefficients<Lexicographical,T,Index> & u,
                     const IndexSet<Index> &LambdaCol);

// Calculates ||u - u_{\Lambda}-f||_H via energy norm ansatz
template <typename T, typename Index, typename MA, typename RHS>
T
computeErrorInH1Norm(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
                     T H1NormOfExactSolution, bool optimized=false);

template <typename T, typename Index, typename SOLVER, typename MA_H, typename RHS_H>
void
postprocessing_H1(SOLVER& Solver, MA_H &A_H1, RHS_H &F_H1, T H1norm, const char* filename);

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_L2(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P);                    

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_H1(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P);
                            
template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_W0T(Coefficients<Lexicographical,T,Index2D> & u,
                            Coefficients<Lexicographical,T,Index2D> & u_exact,
                            const Preconditioner &P);

}    //namespace lawa

#include <lawa/methods/adaptive/postprocessing/error_estimators.tcc>

#endif

