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


#ifndef APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS1D_H
#define APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS1D_H 1

#include <applications/finance/processes/processtypes.h>

namespace lawa {

template < typename T, ProcessType1D Type>
struct ProcessParameters1D
{

};

template <typename T>
struct ProcessParameters1D<T,CGMY>
{
  ProcessParameters1D(T _r, T _k_C, T _k_G, T _k_M, T _k_Y);

  T r;
  T k_C;
  T k_G;
  T k_M;
  T k_Y;
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,CGMY> &processparameters);

template <typename T>
struct ProcessParameters1D<T,CGMYe>
{
  ProcessParameters1D(T _r, T _k_C, T _k_G, T _k_M, T _k_Y, T _sigma);

  T r;
  T k_C;
  T k_G;
  T k_M;
  T k_Y;
  T sigma;
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,CGMYe> &processparameters);

template <typename T>
struct ProcessParameters1D<T,BlackScholes>
{
    ProcessParameters1D(T _r, T _sigma);

    T r;
    T sigma;
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,BlackScholes> &processparameters);

}   // namespace lawa

#include <applications/finance/processes/processparameters1d.tcc>

#endif  // APPLICATIONS_FINANCE_PROCESSES_PROCESSPARAMETERS1D_H
