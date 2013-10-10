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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_THETASCHEMEMULTITREEAWGM_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_THETASCHEMEMULTITREEAWGM_H 1

#include <map>
#include <cstring>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/algorithms/algorithms.h>

namespace lawa {

template <typename Index, typename ThetaTimeStepSolver>
struct ThetaSchemeAWGM {

    typedef typename ThetaTimeStepSolver::T T;

    typedef typename IndexSet<Index>::const_iterator                          const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::iterator          coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator    const_coeff_it;

    ThetaSchemeAWGM(ThetaTimeStepSolver &_multitree_solver);

    void
    setParameters(T _theta, T _timestep, int numOfTimesteps, T _timestep_eps, int _maxiterations,
                  T _init_cgtol, int _strategy);

    void
    applyPreconditioner(Coefficients<Lexicographical,T,Index> &v);

    void
    applyInvPreconditioner(Coefficients<Lexicographical,T,Index> &v);

    void
    solve(Coefficients<Lexicographical,T,Index> &u, int &avDof, int &maxDof, int &terminalDof, int j);

    ThetaTimeStepSolver     &timestep_solver;
    T                       theta, timestep;
    int                     numOfTimesteps;
    T                       timestep_eps;
    int                     maxiterations;
    T                       init_cgtol;
    int                     strategy;
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/thetaschemeawgm.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_THETASCHEMEMULTITREEAWGM_H

