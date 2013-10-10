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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDINDEXSET_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDINDEXSET_H 1

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#else
    #include <ext/hash_set>
#endif

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/bspline.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
struct AlignedIndexSet
{
    #ifdef TRONE
        typedef typename std::tr1::unordered_map<PrincipalIndex,IndexSet<AlignedIndex>,
                                                 index_hashfunction<PrincipalIndex>,
                                                 index_eqfunction<PrincipalIndex> >
                         IndexToIndexSetMap;
    #elif BOOST
        typedef typename boost::unordered_map<PrincipalIndex,IndexSet<AlignedIndex>,
                                                 index_hashfunction<PrincipalIndex>,
                                                 index_eqfunction<PrincipalIndex> >
                         IndexToIndexSetMap;
    #else
        typedef typename __gnu_cxx::hash_map<PrincipalIndex,IndexSet<AlignedIndex>,
                                             index_hashfunction<PrincipalIndex>,
                                             index_eqfunction<PrincipalIndex> >
                         IndexToIndexSetMap;
    #endif

    typedef typename IndexSet<Index>::const_iterator            const_set_index_it;
    typedef typename IndexToIndexSetMap::const_iterator         const_map_prinindex_it;
    typedef typename IndexToIndexSetMap::iterator               map_prinindex_it;

    typedef typename IndexSet<AlignedIndex>::const_iterator     const_set_aligindex_it;
    typedef typename IndexSet<AlignedIndex>::iterator           set_aligindex_it;

    AlignedIndexSet(void);

    AlignedIndexSet(size_t n1, size_t n2);

    /*
     * Aligns the indexset Lambda w.r.t. to the principal index assuming that
     * Index = (PrincipalIndex,AlignedIndex)
     */
    void
    align_x1(const IndexSet<Index> &Lambda);

    /*
     * Aligns the indexset Lambda w.r.t. to the principal index assuming that
     * Index = (AlignedIndex,PrincipalIndex)
     */
    void
    align_x2(const IndexSet<Index> &Lambda);

    void
    clear();

    IndexToIndexSetMap map;

    size_t n1, n2;

};

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
std::ostream& operator<< (std::ostream &s,
                          const AlignedIndexSet<Index,PrincipalIndex,AlignedIndex> &alignedLambda);

}   // namespace

#include <lawa/methods/adaptive/datastructures/alignedindexset.tcc>

#endif  //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ALIGNEDINDEXSET_H
