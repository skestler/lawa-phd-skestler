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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H 1

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#else
    #include <ext/hash_set>
#endif

#include <iostream>
#include <list>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

#define COEFFBYLEVELSIZE 1024

template <typename T>
struct CoefficientsByLevel
{
    #ifdef TRONE
        typedef typename std::tr1::unordered_map<long, T> TranslationIndexToValueMap;
        //typedef typename std::map<long, T> TranslationIndexToValueMap;
    #elif BOOST
        typedef typename boost::unordered_map<long, T> TranslationIndexToValueMap;
    #else
        typedef typename __gnu_cxx::hash_map<long, T> TranslationIndexToValueMap;
    #endif

    typedef typename TranslationIndexToValueMap::const_iterator const_it;
    typedef typename TranslationIndexToValueMap::iterator       iter;
    typedef typename TranslationIndexToValueMap::value_type     val_type;

    CoefficientsByLevel(void);

    CoefficientsByLevel(short _j, size_t n);

    void
    set(short _j, size_t n);

    CoefficientsByLevel<T>&
    operator=(const CoefficientsByLevel<T> &_coeff);

    CoefficientsByLevel<T>&
    operator+=(const CoefficientsByLevel<T> &_coeff);

    void
    setToZero();

    short j;
    TranslationIndexToValueMap map;

};

template <typename T>
std::ostream& operator<<(std::ostream &s, const CoefficientsByLevel<T> &_treecoeff);

template <typename T>
struct TreeCoefficients1D
{
    typedef typename CoefficientsByLevel<T>::val_type                           val_type;
    typedef typename CoefficientsByLevel<T>::const_it                           const_by_level_it;
    typedef typename CoefficientsByLevel<T>::iter                               by_level_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator    const_coeff1d_it;
    typedef typename std::list<const Index1D*>::const_iterator                  const_indexlist_it;

    TreeCoefficients1D(size_t n, int basis_j0);

    TreeCoefficients1D<T>&
    operator=(const TreeCoefficients1D<T> &_coeff);

    TreeCoefficients1D<T>&
    operator=(const Coefficients<Lexicographical,T,Index1D> &_coeff);

    TreeCoefficients1D<T>&
    operator=(const std::list<const Index1D* > &_coeff);

    TreeCoefficients1D<T>&
    operator-=(const Coefficients<Lexicographical,T,Index1D> &_coeff);

    TreeCoefficients1D<T>&
    operator-=(const TreeCoefficients1D<T> &_coeff);

    TreeCoefficients1D<T>&
    operator+=(const Coefficients<Lexicographical,T,Index1D> &_coeff);

    TreeCoefficients1D<T>&
    operator*=(T factor);

    const CoefficientsByLevel<T>&
    operator[](short j)  const;

    CoefficientsByLevel<T>&
    operator[](short j);

    template<typename Index, typename PrincipalIndex, CoordinateDirection CoordX>
    void
    addTo(const PrincipalIndex &lambda, Coefficients<Lexicographical,T,Index> &v);

    void
    setToZero();

    int
    size();

    int
    getMaxTreeLevel();

    int
    setMaxTreeLevel(int j);

    T
    norm(T factor);

    int offset;     // is supposed to be j0-1!!
    CoefficientsByLevel<T> bylevel[JMAX+1];
    int maxTreeLevel;

};

template<typename T>
void
fromTreeCoefficientsToCoefficients(const TreeCoefficients1D<T> &tree_v,
                                 Coefficients<Lexicographical,T,Index1D> &v);

template<typename T>
void
fromCoefficientsToTreeCoefficients(const Coefficients<Lexicographical,T,Index1D> &v,
                                 TreeCoefficients1D<T> &tree_v);

template<typename T, typename ScalingOperator>
void
scale(const TreeCoefficients1D<T> &x, const ScalingOperator &D, TreeCoefficients1D<T> y);

template <typename T>
std::ostream& operator<<(std::ostream &s, const TreeCoefficients1D<T> &_treecoeff);


}   // namespace lawa

#include <lawa/methods/adaptive/datastructures/treecoefficients1d.tcc>

#endif  //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H
