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

namespace lawa {

template <typename Index>
Entry<Index>::Entry(const Index &_index1, const Index &_index2)
: row_index(_index1), col_index(_index2)
{
}

template <typename Index>
std::ostream& operator<<(std::ostream &s, const Entry<Index> &entry) {
    s << "[" << entry.row_index << ", " << entry.col_index  << "]";
    return s;
}

template <typename SortingType>
struct lt<AbsoluteValue, SortingType>
{
    bool operator()(const SortingType &left, const SortingType &right) const
        {
            return (fabs(left) > fabs(right));    //todo: Is this the right call for fabs (template?)
        }
};

} //namespace lawa
