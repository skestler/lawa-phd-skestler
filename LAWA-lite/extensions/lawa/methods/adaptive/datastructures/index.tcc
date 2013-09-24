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

// dim(Index) = dim(Index1) + dim(Index2)
template <>
struct
Split<Index2D,Index1D,Index1D,XOne>{
    inline
    const void operator()(const Index2D &index, Index1D &prinIndex, Index1D &aligIndex) const
    { prinIndex = index.index1; aligIndex = index.index2; return;  }
};

template <>
struct
Join<Index2D,Index1D,Index1D,XOne>{
    inline
    const void operator()(const Index1D &prinIndex, const Index1D &aligIndex, Index2D &index) const
    { index.index1 = prinIndex; index.index2 = aligIndex; return;  }
};

template <>
struct
Split<Index2D,Index1D,Index1D,NotXOne>{
    inline
    const void operator()(const Index2D &index, Index1D &prinIndex, Index1D &aligIndex) const
    { prinIndex = index.index2; aligIndex = index.index1; return;  }
};

template <>
struct
Join<Index2D,Index1D,Index1D,NotXOne>{
    inline
    const void operator()(const Index1D &prinIndex, const Index1D &aligIndex, Index2D &index) const
    { index.index2 = prinIndex; index.index1 = aligIndex; return;  }
};

template <>
struct
Split<Index2D,Index1D,Index1D,XTwo>{
    inline
    const void operator()(const Index2D &index, Index1D &prinIndex, Index1D &aligIndex) const
    { prinIndex = index.index2; aligIndex = index.index1; return;  }
};

template <>
struct
Join<Index2D,Index1D,Index1D,XTwo>{
    inline
    const void operator()(const Index1D &prinIndex, const Index1D &aligIndex, Index2D &index) const
    { index.index2 = prinIndex; index.index1 = aligIndex; return;  }
};

template <>
struct
Split<Index2D,Index1D,Index1D,NotXTwo>{
    inline
    const void operator()(const Index2D &index, Index1D &prinIndex, Index1D &aligIndex) const
    { prinIndex = index.index1; aligIndex = index.index2; return;  }
};

template <>
struct
Join<Index2D,Index1D,Index1D,NotXTwo>{
    inline
    const void operator()(const Index1D &prinIndex, const Index1D &aligIndex, Index2D &index) const
    { index.index1 = prinIndex; index.index2 = aligIndex; return;  }
};

// Three dimensions

template <>
struct
Split<Index3D,Index1D,Index2D,XOne>{
    inline
    const void operator()(const Index3D &index, Index1D &prinIndex, Index2D &aligIndex) const
    { prinIndex = index.index1;
      aligIndex.index1 = index.index2; aligIndex.index2 = index.index3; return;  }
};

template <>
struct
Join<Index3D,Index1D,Index2D,XOne>{
    inline
    const void operator()(const Index1D &prinIndex, const Index2D &aligIndex, Index3D &index) const
    { index.index1 = prinIndex;
      index.index2 = aligIndex.index1; index.index3 = aligIndex.index2; return;  }
};

template <>
struct
Split<Index3D,Index2D,Index1D,NotXOne>{
    inline
    const void operator()(const Index3D &index, Index2D &prinIndex, Index1D &aligIndex) const
    { prinIndex.index1 = index.index2; prinIndex.index2 = index.index3;
      aligIndex = index.index1; return;  }
};

template <>
struct
Join<Index3D,Index2D,Index1D,NotXOne>{
    inline
    const void operator()(const Index2D &prinIndex, const Index1D &aligIndex, Index3D &index) const
    { index.index2 = prinIndex.index1; index.index3 = prinIndex.index2;
      index.index1 = aligIndex; return;  }
};

template <>
struct
Split<Index3D,Index1D,Index2D,XTwo>{
    inline
    const void operator()(const Index3D &index, Index1D &prinIndex, Index2D &aligIndex) const
    { prinIndex = index.index2;
      aligIndex.index1 = index.index1; aligIndex.index2 = index.index3; return;  }
};

template <>
struct
Join<Index3D,Index1D,Index2D,XTwo>{
    inline
    const void operator()(const Index1D &prinIndex, const Index2D &aligIndex, Index3D &index) const
    { index.index2 = prinIndex;
      index.index1 = aligIndex.index1; index.index3 = aligIndex.index2; return;  }
};

template <>
struct
Split<Index3D,Index2D,Index1D,NotXTwo>{
    inline
    const void operator()(const Index3D &index, Index2D &prinIndex, Index1D &aligIndex) const
    { prinIndex.index1 = index.index1; prinIndex.index2 = index.index3;
      aligIndex = index.index2; return;  }
};

template <>
struct
Join<Index3D,Index2D,Index1D,NotXTwo>{
    inline
    const void operator()(const Index2D &prinIndex, const Index1D &aligIndex, Index3D &index) const
    { index.index1 = prinIndex.index1; index.index3 = prinIndex.index2;
      index.index2 = aligIndex; return;  }
};

template <>
struct
Split<Index3D,Index1D,Index2D,XThree>{
    inline
    const void operator()(const Index3D &index, Index1D &prinIndex, Index2D &aligIndex) const
    { prinIndex = index.index3;
      aligIndex.index1 = index.index1; aligIndex.index2 = index.index2; return;  }
};

template <>
struct
Join<Index3D,Index1D,Index2D,XThree>{
    inline
    const void operator()(const Index1D &prinIndex, const Index2D &aligIndex, Index3D &index) const
    { index.index3 = prinIndex;
      index.index1 = aligIndex.index1; index.index2 = aligIndex.index2; return;  }
};

template <>
struct
Split<Index3D,Index2D,Index1D,NotXThree>{
    inline
    const void operator()(const Index3D &index, Index2D &prinIndex, Index1D &aligIndex) const
    { prinIndex.index1 = index.index1; prinIndex.index2 = index.index2;
      aligIndex = index.index3; return;  }
};

template <>
struct
Join<Index3D,Index2D,Index1D,NotXThree>{
    inline
    const void operator()(const Index2D &prinIndex, const Index1D &aligIndex, Index3D &index) const
    { index.index1 = prinIndex.index1; index.index2 = prinIndex.index2;
      index.index3 = aligIndex; return;  }
};


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
