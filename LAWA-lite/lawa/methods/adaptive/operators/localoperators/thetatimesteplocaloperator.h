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
#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_THETATIMESTEPLOCALOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_THETATIMESTEPLOCALOPERATOR_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename Index, typename StiffnessMatrixLocalOperator,
          typename MassMatrixLocalOperator=StiffnessMatrixLocalOperator>
class ThetaTimeStepLocalOperator {

    public:
        typedef typename StiffnessMatrixLocalOperator::T T;
        typedef typename Coefficients<Lexicographical,T,Index>::iterator    coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index>::iterator    const_coeff_it;
        typedef typename IndexSet<Index1D>::const_iterator                  const_set1d_it;

        ThetaTimeStepLocalOperator(T _theta, T _timestep,
                                   StiffnessMatrixLocalOperator &_StiffnessMatrixLocalOp);

        ThetaTimeStepLocalOperator(T _theta, T _timestep,
                                   StiffnessMatrixLocalOperator &_StiffnessMatrixLocalOp,
                                   MassMatrixLocalOperator      &_MassMatrixLocalOp);

        void
        setThetaTimeStepParameters(T _theta, T _timestep);

        void
        evalM(Coefficients<Lexicographical,T,Index> &v,
              Coefficients<Lexicographical,T,Index> &Mv, const char* evalType);

        void
        evalA(Coefficients<Lexicographical,T,Index> &v,
              Coefficients<Lexicographical,T,Index> &Av, const char* evalType);


        void
        eval(const Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av);

        template <typename Preconditioner>
        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, const char* evalType);

        T                               theta;
        T                               timestep;
        bool                            MassIsIdentity;
        StiffnessMatrixLocalOperator    &StiffnessMatrixLocalOp;
        MassMatrixLocalOperator         &MassMatrixLocalOp;
};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/localoperators/thetatimesteplocaloperator.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_THETATIMESTEPLOCALOPERATOR_H
