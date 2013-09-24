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

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {
    
Index1D::Index1D(void)
: j(0), k(0), xtype(XBSpline), val(0), linearindex(0)
{
}

Index1D::Index1D(int _j, int _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype), val(xtype), linearindex(0)
{
    /*
     std::cout << "Index1D.val = " << val << std::endl;
     for (int i=63; i>=0; --i) {
     std::cout << ((val & (1l << i)) ? 1 : 0);
     }
     std::cout<<std::endl;
     val = (val << 16) | j;
     for (int i=63; i>=0; --i) {
     std::cout << ((val & (1l << i)) ? 1 : 0);
     }
     std::cout<<std::endl;
     val = (val << 32 | (unsigned int) k);
     for (int i=63; i>=0; --i) {
     std::cout << ((val & (1l << i)) ? 1 : 0);
     }
     std::cout<<    std::endl;
     */
    val = (((val << 16) | (unsigned short) j) << 32) | (unsigned int) k;
}

Index1D::Index1D(const Index1D &index)
: j(index.j), k(index.k), xtype(index.xtype), val(index.val), linearindex(index.linearindex)
{
}

Index1D::~Index1D(void)
{
}


std::ostream& operator<<(std::ostream &s, const Index1D &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling," << _i.j << "," << _i.k;
    } else {
        s << "wavelet," << _i.j << "," << _i.k;
    }
    return s;
}


Index2D::Index2D(const Index1D &_index1, const Index1D &_index2)
: linearindex(0), index1(_index1), index2(_index2)
{
}

Index2D::~Index2D(void)
{
}

std::ostream& operator<<(std::ostream &s, const Index2D &_i)
{
    s << _i.index1 << "," << _i.index2;
    return s;
}


Index3D::Index3D(const Index1D &_index1, const Index1D &_index2, const Index1D &_index3)
: linearindex(0), index1(_index1), index2(_index2), index3(_index3)
{
}

Index3D::~Index3D(void)
{
}

std::ostream& operator<<(std::ostream &s, const Index3D &_i)
{
    s <<  _i.index1 << "," << _i.index2 << "," << _i.index3;
    return s;
}

//Bitmask implementation
bool 
lt<Lexicographical, Index1D>::operator()(const Index1D &left, const Index1D &right) const
{
    return left.val < right.val;
}

bool
lt<Lexicographical, Index1D>::operator()(const Entry<Index1D> &left, const Entry<Index1D> &right) const
{
    // sort Operator row-wise
    if (left.row_index.val != right.row_index.val) return left.row_index.val < right.row_index.val;
    else                                           return left.col_index.val < right.col_index.val;
}

bool
lt<Lexicographical, Index2D>::operator()(const Index2D &left, const Index2D &right) const
{
    if (left.index1.val != right.index1.val) return left.index1.val < right.index1.val;
    else                                     return left.index2.val < right.index2.val;
}
    
bool 
lt<Lexicographical, Index2D>::operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const
{
    // sort Operator row-wise
    if (left.row_index.index1.val != right.row_index.index1.val) {
        return left.row_index.index1.val < right.row_index.index1.val;
    }
    else if (left.row_index.index2.val != right.row_index.index2.val) {
        return left.row_index.index2.val < right.row_index.index2.val;
    }
    if (left.col_index.index1.val != right.col_index.index1.val) {
        return left.col_index.index1.val < right.col_index.index1.val;
    }
    else {
        return left.col_index.index2.val < right.col_index.index2.val;
    }
}
    
bool
lt<Lexicographical, Index3D>::operator()(const Index3D &left, const Index3D &right) const
{
    if (left.index1.val != right.index1.val)         return left.index1.val < right.index1.val;
    else if (left.index2.val != right.index2.val)   return left.index2.val < right.index2.val;
    else                                            return left.index3.val < right.index3.val;
}
    
bool
lt<Lexicographical, Index3D>::operator()(const Entry<Index3D> &left, const Entry<Index3D> &right) const
{
    // sort Operator row-wise
    if (left.row_index.index1.val != right.row_index.index1.val) {
        return left.row_index.index1.val < right.row_index.index1.val;
    }
    else if (left.row_index.index2.val != right.row_index.index2.val) {
        return left.row_index.index2.val < right.row_index.index2.val;
    }
    else if (left.row_index.index3.val != right.row_index.index3.val) {
        return left.row_index.index3.val < right.row_index.index3.val;
    }
    if (left.col_index.index1.val != right.col_index.index1.val) {
        return left.col_index.index1.val < right.col_index.index1.val;
    }
    else if (left.col_index.index2.val != right.col_index.index2.val) {
        return left.col_index.index2.val < right.col_index.index2.val;
    }
    else {
        return left.col_index.index3.val < right.col_index.index3.val;
    }
}

bool
index_eqfunction<Index1D>::operator()(const Index1D& leftindex, const Index1D& rightindex) const
{
    return (leftindex.val == rightindex.val);
}

bool
index_eqfunction<Index2D>::operator()(const Index2D& leftindex, const Index2D& rightindex) const
{
    return (    (leftindex.index1.val == rightindex.index1.val)
            && (leftindex.index2.val == rightindex.index2.val)  );
}

size_t
index_hashfunction<Index1D>::operator()(const Index1D& index) const
{
    boost::hash<long> hasher;
    return hasher(index.val);
//    return index.val % 26681;
}

size_t
index_hashfunction<Index2D>::operator()(const Index2D& index) const
{
    std::size_t seed = 0;
    boost::hash_combine(seed, index.index1.val);
    boost::hash_combine(seed, index.index2.val);
    return seed;
//    return (3*index.index1.val + 7*index.index2.val) % 999331;
}

bool
entry_eqfunction<Index1D>::operator()(const Entry<Index1D>& leftentry,
                                      const Entry<Index1D>& rightentry) const
{
    return (    (leftentry.col_index.val == rightentry.col_index.val)
             && (leftentry.row_index.val == rightentry.row_index.val) );
}

size_t
entry_hashfunction<Index1D>::operator()(const Entry<Index1D>& entry) const
{
    std::size_t seed = 0;
    boost::hash_combine(seed, entry.col_index.val);
    boost::hash_combine(seed, entry.row_index.val);
    return seed;
}


} //namespace lawa
