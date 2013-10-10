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

#ifndef APPLICATIONS_FINANCE_FOURIERPRICING_FOURIERPRICER1D_H
#define APPLICATIONS_FINANCE_FOURIERPRICING_FOURIERPRICER1D_H 1


#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_interp.h>

#include <applications/finance/processes/processtypes.h>
#include <applications/finance/fourierpricing/cgmycharacteristicfunction1d.h>
#include <applications/finance/fourierpricing/characteristicfunction1d.h>

namespace lawa {

template <typename T, ProcessType1D PType>
struct FourierPricer1D
{

    FourierPricer1D(CharacteristicFunction1D<T,PType> &_charfunc,
                    T _S0, T _maturity, T _K1, T _K2);

    void solve(T A, int J); //integral bound A, 2^J function evaluations for the FFT

    T
    operator()(T K);

    gsl_complex
    zeta(T v);

    CharacteristicFunction1D<T,PType>      &charfunc;
    T                                      S0, maturity, r;  //spot S0, interest rate r, maturity
    T                                      K1, K2;           //range of strike prices [K1,K2]
    gsl_interp                             *interpol;
    gsl_interp_accel                       *accel;
    double                                 *xa,*ya;

};

}   // namespace lawa

#include <applications/finance/fourierpricing/fourierpricer1d.tcc>

#endif  // APPLICATIONS_FINANCE_FOURIERPRICING_FOURIERPRICER1D_H
