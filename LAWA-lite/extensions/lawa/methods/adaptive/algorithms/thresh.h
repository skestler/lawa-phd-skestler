#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    ABSOLUTE_THRESH(const Coefficients<Lexicographical,T,Index > &v, T eta);

template <typename T>
    Coefficients<Lexicographical,T,Index1D >
    THRESH(const Coefficients<Lexicographical,T,Index1D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T>
    Coefficients<Lexicographical,T,Index2D >
    THRESH(const Coefficients<Lexicographical,T,Index2D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T>
    Coefficients<Lexicographical,T,Index2D >
    MULTITREE_THRESH(const Coefficients<Lexicographical,T,Index2D > &v, T eta);

template <typename T>
    Coefficients<Lexicographical,T,Index3D >
    THRESH(const Coefficients<Lexicographical,T,Index3D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    THRESH(const Coefficients<AbsoluteValue,T,Index > &v, T eta, bool hp=false);



} // namespace lawa

#include <lawa/methods/adaptive/algorithms/thresh.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H

