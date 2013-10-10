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

#ifndef LAWA_SETTINGS_TYPETRAITS_H
#define LAWA_SETTINGS_TYPETRAITS_H 1

#include <lawa/settings/enum.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

//--- IsPrimal

template <typename X>
struct IsPrimal
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
          template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsPrimal<Some<T, Primal, Domain, Cons> >
{
    static const bool value = true;
};

//--- IsOrthogonal

template <typename X>
struct IsOrthogonal
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
    template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsOrthogonal<Some<T, Orthogonal, Domain, Cons> >
{
    static const bool value = true;
};

//--- IsDual

template <typename X>
struct IsDual
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
          template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsDual<Some<T, Dual, Domain, Cons> >
{
    static const bool value = true;
};

//--- PrimalOrDual
template <typename X>
struct PrimalOrDual
{
    static const bool value = (IsPrimal<X>::value || IsDual<X>::value);
};

//--- PrimalOrDualOrOrthogonal
template <typename X>
struct PrimalOrDualOrOrthogonal
{
    static const bool value = (IsPrimal<X>::value || IsDual<X>::value || IsOrthogonal<X>::value);
};

//--- BothPrimal
template <typename X, typename Y>
struct BothPrimal
{
    static const bool value = IsPrimal<X>::value && IsPrimal<Y>::value;
};

//--- BothDual
template <typename X, typename Y>
struct BothDual
{
    static const bool value = IsDual<X>::value && IsDual<Y>::value;
};

//--- BothOrthogonal
template <typename X, typename Y>
struct BothOrthogonal
{
    static const bool value = IsOrthogonal<X>::value && IsOrthogonal<Y>::value;
};

//--- PrimalOrDual
template <typename X, typename Y>
struct PrimalAndDual
{
    static const bool value = (IsPrimal<X>::value && IsDual<Y>::value);
};

//--- IsWavelet
template <typename X>
struct IsWavelet
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct IsWavelet<Wavelet<T,Side,Domain,Cons> >
{
    static const bool value = true;
};

//--- IsBSpline
template <typename X>
struct IsBSpline
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct IsBSpline<BSpline<T,Side,Domain,Cons> >
{
    static const bool value = true;
};

//--- IsPeriodic
template <typename X>
struct IsPeriodic
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsPeriodic<Wavelet<T,Side,Periodic,Cons> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsPeriodic<BSpline<T,Side,Periodic,Cons> >
{
    static const bool value = true;
};

//--- IsRealline
template <typename X>
struct IsRealline
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsRealline<Wavelet<T,Side,R,Cons> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsRealline<BSpline<T,Side,R,Cons> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsRealline<Basis<T,Side,R,Cons> >
{
    static const bool value = true;
};

//--- IsSparseMulti
template <typename X>
struct IsSparseMulti
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, DomainType Domain>
struct IsSparseMulti<Wavelet<T,Side,Domain,SparseMulti> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, DomainType Domain>
struct IsSparseMulti<BSpline<T,Side,Domain,SparseMulti> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, DomainType Domain>
struct IsSparseMulti<Basis<T,Side,Domain,SparseMulti> >
{
    static const bool value = true;
};

//--- IsIndex1D
template <typename X>
struct IsIndex1D
{
    static const bool value = false;
};

template <>
struct IsIndex1D<Index1D>
{
    static const bool value = true;
};

//--- IsIndex2D
template <typename X>
struct IsIndex2D
{
    static const bool value = false;
};

template <>
struct IsIndex2D<Index2D>
{
    static const bool value = true;
};

//--- IsIndex3D
template <typename X>
struct IsIndex3D
{
    static const bool value = false;
};

template <>
struct IsIndex3D<Index3D>
{
    static const bool value = true;
};

} // namespace lawa
 
#endif // LAWA_SETTINGS_TYPETRAITS_H
