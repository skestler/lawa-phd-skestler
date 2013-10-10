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
#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H 1

#include <iostream>
#include <vector>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/algorithms/linearsystemsolvers.h>
#include <lawa/methods/adaptive/postprocessing/postprocessing.h>

namespace lawa {

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
class S_ADWAV {
    public:
        S_ADWAV(const Basis &basis, MA &A, RHS &F, T contraction, T start_threshTol,
                T _linTol=1e-6, T _resTol=1e-4, int _NumOfIterations=10, int _MaxItsPerThreshTol=5,
                T _eps=1e-2, int _MaxSizeLambda = 400, T _resStopTol=0.1,
                std::vector<int> _Jmaxvec = std::vector<int>(0));

        //solver for symmetric elliptic problems
        void solve(const IndexSet<Index> &Initial_Lambda, const char *linsolvertype,
                   const char *filename, int assemble_matrix=2, T H1norm=0.);
        //solver for symmetric elliptic problems
        void solve_cg(const IndexSet<Index> &Initial_Lambda, int assemble_matrix=1, T H1norm=0.);
        //solver for symmetric elliptic problems without B-Splines
        void solve_cg_WO_XBSpline(const IndexSet<Index> &Initial_Lambda, int assemble_matrix=1, T H1norm=0.);
        //solver for elliptic problems
        void solve_gmres(const IndexSet<Index> &Initial_Lambda, int assemble_matrix=1);
        void solve_gmresm(const IndexSet<Index> &Initial_Lambda, int assemble_matrix=1);
        //solver for indefinite problems
        void solve_cgls(const IndexSet<Index> &Initial_Lambda, int assemble_matrix=1);
        
        void
        set_parameters(T _contraction, T _threshTol, T _linTol=1e-6, T _resTol=1e-4, 
                       int _NumOfIterations=10, int _MaxItsPerThreshTol=5, T _eps=1e-2, 
                       int _MaxSizeLambda = 400, T _resStopTol=0.1, 
                       std::vector<int> _Jmaxvec = std::vector<int>(0));
        void
        get_parameters(T& _contraction, T& _threshTol, T& _linTol, T& _resTol, 
                       int& _NumOfIterations, int& _MaxItsPerThreshTol, T& _eps, 
                       int& _MaxSizeLambda, T& _resStopTol,
                       std::vector<int>& _Jmaxvec);
    
        std::vector<Coefficients<Lexicographical,T,Index> > solutions;
        std::vector<T>               residuals;
        std::vector<T>               times;
        std::vector<T>               linsolve_iterations;
        std::vector<T>               toliters;

        const Basis &basis;
    
    private:
        MA &A;
        RHS &F;
        T contraction, threshTol, linTol, resTol;
        int NumOfIterations; 
        int MaxItsPerThreshTol;
        T eps;
        int MaxSizeLambda;
        T resStopTol;
        std::vector<int> Jmaxvec;
        

};

} //namespace lawa

#include <lawa/methods/adaptive/solvers/s_adwav.tcc>

#endif //LAWA_METHODS_ADAPTIVE_SOLVERS_S_ADWAV_H

