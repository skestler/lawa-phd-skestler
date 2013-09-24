/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_AUX_INTEGRAL_VS_GMP_H
#define LAWA_AUX_INTEGRAL_VS_GMP_H 1


#include <complex>
//#include <gmpxx.h>
//#include <extensions/gmpmath/mpreal.h>

namespace lawa {

//------------------------------------------------------------------------------

template <typename T>
struct IntegralReal
{
};

template <>
struct IntegralReal<float>
{
    typedef float Type;
};

template <>
struct IntegralReal<double>
{
    typedef double Type;
};

template <typename T>
struct IntegralInt
{
};

template <>
struct IntegralInt<int>
{
    typedef int Type;
};

template <>
struct IntegralInt<long>
{
    typedef long Type;
};

//------------------------------------------------------------------------------

/*
template <typename T>
struct GMPReal {
};

template <>
struct GMPReal<mpfr::mpreal>
{
    typedef mpfr::mpreal Type;
};

template <typename T>
struct GMPInt {
};

template <>
struct GMPInt<mpz_class>
{
    typedef mpz_class Type;
};
*/

} // namespace lawa

#endif // LAWA_AUX_INTEGRAL_VS_GMP_H


