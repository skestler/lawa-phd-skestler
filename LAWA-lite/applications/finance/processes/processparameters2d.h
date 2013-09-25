/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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


#ifndef APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS2D_H
#define APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS2D_H 1

#include <applications/finance/processes/processtypes.h>

namespace lawa {

template < typename T, ProcessType2D Type>
struct ProcessParameters2D
{

};

template <typename T>
struct ProcessParameters2D<T,BlackScholes2D>
{
  ProcessParameters2D(T _r, T _sigma1, T _sigma2, T _rho, T _u11, T _u12, T _u21, T _u22);

  T r;
  T sigma1, sigma2;
  T rho;
  T u11, u12, u21, u22;
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters2D<T,BlackScholes2D> &processparameters);

template <typename T>
struct ProcessParameters2D<T,CGMYeUnivariateJump2D>
{
  ProcessParameters2D(T _r, T _sigma1, T _sigma2, T _rho, T _k_C1, T _k_G1, T _k_M1, T _k_Y1,
                      T _k_C2, T _k_G2, T _k_M2, T _k_Y2);

  T r;
  T sigma1, sigma2;
  T rho;
  T k_C1, k_G1, k_M1, k_Y1;
  T k_C2, k_G2, k_M2, k_Y2;
  ProcessParameters1D<T,CGMYe> proc_param1, proc_param2;
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters);


}   // namespace lawa

#include <applications/finance/processes/processparameters2d.tcc>

#endif  // APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS2D_H
