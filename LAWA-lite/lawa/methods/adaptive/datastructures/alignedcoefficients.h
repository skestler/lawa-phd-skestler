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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDCOEFFICIENTS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDCOEFFICIENTS_H 1

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#else
    #include <ext/hash_map>
#endif

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection X>
struct AlignedCoefficients
{
    #ifdef TRONE
        typedef typename std::tr1::unordered_map<PrincipalIndex,
                                                 Coefficients<Lexicographical,T,AlignedIndex>,
                                                 index_hashfunction<PrincipalIndex>,
                                                 index_eqfunction<PrincipalIndex> >
                         IndexToCoefficientsMap;
    #elif BOOST
        typedef typename boost::unordered_map<PrincipalIndex,
                                                 Coefficients<Lexicographical,T,AlignedIndex>,
                                                 index_hashfunction<PrincipalIndex>,
                                                 index_eqfunction<PrincipalIndex> >
                         IndexToCoefficientsMap;
    #else
        typedef typename __gnu_cxx::hash_map<PrincipalIndex,
                                             Coefficients<Lexicographical,T,AlignedIndex>,
                                             index_hashfunction<PrincipalIndex>,
                                             index_eqfunction<PrincipalIndex> >
                         IndexToCoefficientsMap;
    #endif

    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator        const_coeff_index_it;
    typedef typename IndexToCoefficientsMap::const_iterator                       const_map_prindex_it;
    typedef typename IndexToCoefficientsMap::iterator                             map_prinindex_it;

    typedef typename Coefficients<Lexicographical,T,AlignedIndex>::const_iterator const_coeff_aligindex_it;
    typedef typename Coefficients<Lexicographical,T,AlignedIndex>::iterator       coeff_aligindex_it;

    AlignedCoefficients(void);

    AlignedCoefficients(size_t n1, size_t n2);

    /// Case 1: PrincipalIndex = Index1D, AlignedIndex = Index"(n-1)"D -> Index = Index"n"D
    ///         Let now Index = (\lambda_1,...,\lambda_n). Then the principal index and the
    ///         aligned index are determined by CoordX. E.g. when CoordX=XTwo, we get
    ///         principal index = \lambda_2, aligned index = (\lambda_1,\lambda_3,...,\lambda_n)
    /// Case 2: PrincipalIndex = Index"n-1"D, AlignedIndex = Index1D -> Index = Index"n"D
    ///         Let now Index = (\lambda_1,...,\lambda_n). Then the principal index and the
    ///         aligned index are determined by CoordX. E.g. when CoordX=XTwo, we get
    ///         principal index = (\lambda_1,\lambda_3,...,\lambda_n), aligned index = \lambda_2

    void
    align(const Coefficients<Lexicographical,T,Index> &coeff, short J=0);

    void
    align(const Coefficients<Lexicographical,T,Index> &coeff, IndexSet<PrincipalIndex> &prinIndices,
          short J=0);

    void
    align_ExcludeAndOrthogonal(const Coefficients<Lexicographical,T,Index> &coeff,
                               const Coefficients<Lexicographical,T,Index> &exclude,
                               const IndexSet<PrincipalIndex> &orthogonal, short J=0);

    IndexToCoefficientsMap map;

    size_t n1, n2;

};

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
struct AlignedCoefficients2
{
    typedef Coefficients<Lexicographical,int,PrincipalIndex>    PrincipalIndexMap;
    typedef typename std::list<const AlignedIndex* >            AlignedIndices;
    typedef typename std::vector<AlignedIndices>                PrincipalIndexToAlignedIndices;

    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator   const_coeff_index_it;
    typedef typename PrincipalIndexMap::const_iterator                       const_coeff_prinindex_it;


    AlignedCoefficients2(void);

    void
    align_x1(const Coefficients<Lexicographical,T,Index> &coeff);

    PrincipalIndexMap               principalIndices;
    PrincipalIndexToAlignedIndices  principalIndexToAlignedIndices;

};

}   // namespace lawa

#include <lawa/methods/adaptive/datastructures/alignedcoefficients.tcc>

#endif  //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDCOEFFICIENTS_H
