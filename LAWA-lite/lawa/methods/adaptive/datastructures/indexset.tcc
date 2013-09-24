/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#include <lawa/math/math.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename Index>
IndexSet<Index>::IndexSet(void)
{
}

#ifdef TRONE
template <typename Index>
IndexSet<Index>::IndexSet(size_t n)
: std::tr1::unordered_set<Index, index_hashfunction<Index>, index_eqfunction<Index> >::unordered_set(n)
{
}
#elif BOOST
template <typename Index>
IndexSet<Index>::IndexSet(size_t n)
: boost::unordered_set<Index, index_hashfunction<Index>, index_eqfunction<Index> >::unordered_set(n)
{
}
#else
template <typename Index>
IndexSet<Index>::IndexSet(size_t n)
:__gnu_cxx::hash_set<Index, index_hashfunction<Index>, index_eqfunction<Index> >::hash_set(n)
{
}
#endif

template <typename Index>
IndexSet<Index>&
IndexSet<Index>::operator=(const IndexSet<Index> &_set)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    this->erase(IndexSet<Index>::begin(), IndexSet<Index>::end());
    if (_set.size() > 0) {
        for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
            this->insert(*lambda);
        }
    }
    return *this;
}

template <typename Index>
IndexSet<Index>
IndexSet<Index>::operator+(const IndexSet<Index> &_set) const
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    IndexSet<Index> ret = *this;
    for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
        ret.insert(*lambda);
    }
    return ret;
}

template <typename Index>
IndexSet<Index> &
IndexSet<Index>::operator+=(const IndexSet<Index> &_set)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
        (*this).insert(*lambda);
    }
    return (*this);
}

template <typename Index>
std::ostream& operator<< (std::ostream &s, const IndexSet<Index> &i)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    s << std::endl << "IndexSet<Index>:" << std::endl;
    if (i.size() > 0) {
        for (const_it lambda = i.begin(); lambda != i.end(); ++lambda) {
            s << "  [" << (*lambda) << "], " << std::endl;
        }
    }
    return s << std::endl;
}

} // namespace lawa

