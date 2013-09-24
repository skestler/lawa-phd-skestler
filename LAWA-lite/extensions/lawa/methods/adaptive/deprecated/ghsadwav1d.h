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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV1D_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV1D_H 1

#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/algorithms/apply1d.h>

namespace lawa {

template <typename T, typename Basis, typename APPLY1D, typename RHS>
class GHS_ADWAV1D {

        typedef typename IndexSet<Index1D>::iterator                             set_it;
        typedef typename IndexSet<Index1D>::const_iterator                       const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff_it;
        typedef typename Coefficients<AbsoluteValue,T,Index1D >::const_iterator  const_coeff_abs_it;
        typedef typename Coefficients<Lexicographical,T,Index1D>::value_type     val_type;

        typedef typename APPLY1D::MAType MA;


    public:

        GHS_ADWAV1D(const Basis &_basis, APPLY1D &_Apply, RHS &_F);

        Coefficients<Lexicographical,T,Index1D>
        SOLVE(T nuM1, T _eps, int NumOfIterations=100, T H1norm=0.);

        std::vector<Coefficients<Lexicographical,T,Index1D> > solutions;
        std::vector<T>                                        residuals;
        std::vector<T>                                        times;
        std::vector<int>                                      linsolve_iterations;


    private:
        const Basis &basis;
        APPLY1D &Apply;
        MA &A;
        RHS &F;
        T cA, CA, kappa;
        T alpha, omega, gamma, theta;
        T eps;

        IndexSet<Index1D>
        GROW(const Coefficients<Lexicographical,T,Index1D> &w, T nu_bar, T &nu);

        Coefficients<Lexicographical,T,Index1D>
        GALSOLVE(const IndexSet<Index1D> &Lambda, const Coefficients<Lexicographical,T,Index1D> &g,
                 const Coefficients<Lexicographical,T,Index1D> &w, T delta, T tol);
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/ghsadwav1d.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV1D_H

