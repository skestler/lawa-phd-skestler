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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H 1

#include <boost/functional/hash.hpp>

#include <lawa/settings/enum.h>
#include <iostream>

namespace lawa {

static boost::hash<long int> hash_long;

//#define JMINOFFSET                 6
#define JMINOFFSET                  6
#define JMAX                       40
#define SIZEHASHINDEX1D         12869//196613
#define SIZELARGEHASHINDEX1D    72869//1572869
#define SIZEHASHINDEX2D       6291469
#define SIZELARGEHASHINDEX2D  6291469
#define SIZELARGEHASHINDEX3D  6291469


struct Index1D
{
    short j;
    long k;
    XType xtype;

    Index1D(void);
    Index1D(int j, long k, XType _xtype);
    Index1D(const Index1D &index);
    ~Index1D();
    short levelSum() const;
};

std::ostream& operator<<(std::ostream &s, const Index1D &_Index);

struct Index2D
{
    Index2D(void);
    Index2D(const Index1D &index1, const Index1D &index2);
    ~Index2D();
    short levelSum() const;

    Index1D index1, index2;

};

std::ostream& operator<<(std::ostream &s, const Index2D &_Index);

struct Index3D
{
    Index3D(void);
    Index3D(const Index1D &index1, const Index1D &index2, const Index1D &index3);
    ~Index3D();
    short levelSum() const;

    Index1D index1, index2, index3;

};

std::ostream& operator<<(std::ostream &s, const Index3D &_Index);

template <typename Index, typename PrincipalIndex, typename AlignedIndex, CoordinateDirection CoordX>
struct
Split{ };

template <typename Index, typename PrincipalIndex, typename AlignedIndex, CoordinateDirection CoordX>
struct
Join{ };


template <typename Index>
struct Entry
{
    Entry(const Index &row_index, const Index &col_index);
    const Index row_index, col_index;    //todo: no copy, but only a reference possible ?!
};

template <typename Index>
std::ostream& operator<<(std::ostream &s, const Entry<Index> &entry);

template <SortingCriterion S, typename SortingType>
struct lt
{
};

/*
 * Compare (lt) operators.
 */

template<>
struct lt<Lexicographical, Index1D>
{
    inline
    bool operator()(const Index1D &leftindex, const Index1D &rightindex) const
    {
        if (leftindex.j!=rightindex.j) {
            return (leftindex.j<rightindex.j);
        } else {
            if ((leftindex.xtype==XBSpline)&&(rightindex.xtype==XWavelet)) {
                return true;
            } else if ((leftindex.xtype==XWavelet)&&(rightindex.xtype==XBSpline)) {
                return false;
            } else {
                return (leftindex.k<rightindex.k);
            }
        }
    }

    bool operator()(const Entry<Index1D> &left, const Entry<Index1D> &right) const;
};

template <>
struct lt<Lexicographical, Index2D>
{
    inline
    bool operator()(const Index2D &left, const Index2D &right) const
    {
        if         ((left.index1.xtype==XBSpline)&&(right.index1.xtype==XWavelet)) {
            return true;
        }
        else if ((left.index1.xtype==XWavelet)&&(right.index1.xtype==XBSpline)) {
            return false;
        }
        else if ((left.index2.xtype==XBSpline)&&(right.index2.xtype==XWavelet)) {
            return true;
        }
        else if ((left.index2.xtype==XWavelet)&&(right.index2.xtype==XBSpline)) {
            return false;
        }
        else if (left.index1.j!=right.index1.j) {
            return (left.index1.j<right.index1.j);
        }
        else if (left.index1.k!=right.index1.k) {
            return (left.index1.k<right.index1.k);
        }
        else if (left.index2.j!=right.index2.j) {
            return (left.index2.j<right.index2.j);
        }
        else {
            return (left.index2.k<right.index2.k);
        }
    }

    //bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const;
};


/*
 * Equal functions.
 */

template <typename Index>
struct index_eqfunction
{
};

template <>
struct index_eqfunction<Index1D>
{
    inline
    bool operator()(const Index1D& leftindex, const Index1D& rightindex) const
    {
        if (leftindex.k != rightindex.k) return false;
        else {
            int leftval = leftindex.xtype;
            leftval = (((leftval << 16) | (unsigned short) leftindex.j));
            int rightval = rightindex.xtype;
            rightval = (((rightval << 16) | (unsigned short) rightindex.j));
            return (leftval==rightval);
        }
        /*
        if (leftindex.k != rightindex.k) return false;
        else {
            if (leftindex.j != rightindex.j) return false;
            else                             return (leftindex.xtype == rightindex.xtype);
        }
        */
    }
};

template <>
struct index_eqfunction<Index2D>
{
    inline
    bool operator()(const Index2D& leftindex, const Index2D& rightindex) const
    {
        if (leftindex.index1.k != rightindex.index1.k) return false;
        else {
            if (leftindex.index2.k != rightindex.index2.k) return false;
            else {
                int leftval1 = leftindex.index1.xtype;
                leftval1 = (((leftval1 << 16) | (unsigned short) leftindex.index1.j));
                int rightval1 = rightindex.index1.xtype;
                rightval1 = (((rightval1 << 16) | (unsigned short) rightindex.index1.j));
                if (leftval1 != rightval1) return false;
                else {
                    int leftval2 = leftindex.index2.xtype;
                    leftval2 = (((leftval2 << 16) | (unsigned short) leftindex.index2.j));
                    int rightval2 = rightindex.index2.xtype;
                    rightval2 = (((rightval2 << 16) | (unsigned short) rightindex.index2.j));
                    return (leftval2==rightval2);
                }
            }
        }
        //return (leftindex.index1.k == rightindex.index1.k && leftindex.index2.k == rightindex.index2.k &&
        //        leftindex.index1.j == rightindex.index1.j && leftindex.index2.j == rightindex.index2.j &&
        //        leftindex.index1.xtype == rightindex.index1.xtype && leftindex.index2.xtype == rightindex.index2.xtype);

    }
};

template <>
struct index_eqfunction<Index3D>
{
    inline
    bool operator()(const Index3D& leftindex, const Index3D& rightindex) const
    {
        if (leftindex.index1.k != rightindex.index1.k) return false;
        if (leftindex.index2.k != rightindex.index2.k) return false;
        if (leftindex.index3.k != rightindex.index3.k) return false;

        int leftval1 = leftindex.index1.xtype;
        leftval1 = (((leftval1 << 16) | (unsigned short) leftindex.index1.j));
        int rightval1 = rightindex.index1.xtype;
        rightval1 = (((rightval1 << 16) | (unsigned short) rightindex.index1.j));
        if (leftval1 != rightval1) return false;

        int leftval2 = leftindex.index2.xtype;
        leftval2 = (((leftval2 << 16) | (unsigned short) leftindex.index2.j));
        int rightval2 = rightindex.index2.xtype;
        rightval2 = (((rightval2 << 16) | (unsigned short) rightindex.index2.j));
        if (leftval2 != rightval2) return false;

        int leftval3 = leftindex.index3.xtype;
        leftval3 = (((leftval3 << 16) | (unsigned short) leftindex.index3.j));
        int rightval3 = rightindex.index3.xtype;
        rightval3 = (((rightval3 << 16) | (unsigned short) rightindex.index3.j));
        return (leftval3==rightval3);
    }

    //return (leftindex.index1.k == rightindex.index1.k && leftindex.index2.k == rightindex.index2.k && leftindex.index3.k == rightindex.index3.k &&
    //        leftindex.index1.j == rightindex.index1.j && leftindex.index2.j == rightindex.index2.j && leftindex.index3.j == rightindex.index3.j &&
    //        leftindex.index1.xtype == rightindex.index1.xtype && leftindex.index2.xtype == rightindex.index2.xtype  && leftindex.index3.xtype == rightindex.index3.xtype);
};

template <typename Index>
struct entry_eqfunction
{
};

template <>
struct entry_eqfunction<Index1D>
{
    inline
    bool operator()(const Entry<Index1D>& leftentry, const Entry<Index1D>& rightentry) const
    {
        if (leftentry.col_index.k!=rightentry.col_index.k) return false;
        else {
            if (leftentry.row_index.k!=rightentry.row_index.k) return false;
            else {
                if (leftentry.col_index.j!=rightentry.col_index.j) return false;
                else {
                    if (leftentry.row_index.j!=rightentry.row_index.j) return false;
                    else {
                        if (leftentry.col_index.xtype!=rightentry.col_index.xtype) return false;
                        else return (leftentry.row_index.xtype==rightentry.row_index.xtype);
                    }
                }
            }
        }

        /*
        return (leftentry.col_index.k     == rightentry.col_index.k     && leftentry.row_index.k     == rightentry.row_index.k &&
                leftentry.col_index.j     == rightentry.col_index.j     && leftentry.row_index.j     == rightentry.row_index.j &&
                leftentry.col_index.xtype == rightentry.col_index.xtype && leftentry.row_index.xtype == rightentry.row_index.xtype);
        */
    }
};


/*
 * Hash functions.
 */

template <typename Index>
struct index_hashfunction
{
};

template <>
struct index_hashfunction<Index1D>
{
    inline
    size_t operator()(const Index1D& index) const
    {
        // Note: hash_values is taken mod "length of hashtable" automatically!!
        /*
        long pow2ij = (1L << (index.j+JMINOFFSET+index.xtype));
        size_t hash_value = (pow2ij + index.k);

        return hash_value;
        */

        int val = index.xtype;
        val = (((val << 16) | (unsigned short) index.j));

        std::size_t hash_value = 0;
        boost::hash_combine(hash_value, val);
        boost::hash_combine(hash_value, index.k);
        return hash_value;

    }
};

template <>
struct index_hashfunction<Index2D>
{

    // performs better without storing 2^l values... why??
    index_hashfunction(void)
    {
        for (int i=0; i<JMAX+JMINOFFSET+2; ++i) {
            power2i[i] = 1L << i;
        }
    }
    size_t power2i[JMAX+JMINOFFSET+2];


    inline
    size_t operator()(const Index2D& index) const
    {
        //size_t l1 = (1L << (index.index1.j+index.index1.xtype+JMINOFFSET) ) + index.index1.k;
        //size_t l2 = (1L << (index.index2.j+index.index2.xtype+JMINOFFSET) ) + index.index2.k;
        size_t l1 = power2i[index.index1.j+index.index1.xtype+JMINOFFSET] + index.index1.k;
        size_t l2 = power2i[index.index2.j+index.index2.xtype+JMINOFFSET] + index.index2.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t P=SIZELARGEHASHINDEX2D, twoP=2*SIZELARGEHASHINDEX2D;
        return (((((s2+1)%(twoP)) * (s2 % twoP)) % twoP)/2 + s1 % P) % P;
    }

    /*
    inline
    size_t operator()(const Index2D& index) const
    {
        int val1 = index.index1.xtype;
        int val2 = index.index2.xtype;
        val1 = (((val1 << 16) | (unsigned short) index.index1.j));
        val2 = (((val2 << 16) | (unsigned short) index.index2.j));

        std::size_t hash_value = 0;
        boost::hash_combine(hash_value, val1);
        boost::hash_combine(hash_value, val2);
        boost::hash_combine(hash_value, index.index1.k);
        boost::hash_combine(hash_value, index.index2.k);

        return hash_value;
    }
    */
};

template <>
struct index_hashfunction<Index3D>
{

    // performs better without storing 2^l values... why??
    index_hashfunction(void)
    {
        for (int i=0; i<JMAX+JMINOFFSET+2; ++i) {
            power2i[i] = 1L << i;
        }
    }
    size_t power2i[JMAX+JMINOFFSET+2];


    inline
    size_t operator()(const Index3D& index) const
    {
        //size_t l1 = (1L << (index.index1.j+index.index1.xtype+JMINOFFSET) ) + index.index1.k;
        //size_t l2 = (1L << (index.index2.j+index.index2.xtype+JMINOFFSET) ) + index.index2.k;
        size_t l1 = power2i[index.index1.j+index.index1.xtype] + index.index1.k;
        size_t l2 = power2i[index.index2.j+index.index2.xtype] + index.index2.k;
        size_t l3 = power2i[index.index3.j+index.index3.xtype] + index.index3.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t s3 = l1+l2+l3;
        size_t P=SIZELARGEHASHINDEX2D, TwoP=2*SIZELARGEHASHINDEX2D, SixP=6*SIZELARGEHASHINDEX2D;
        return (   ( ( ((s3+2)%SixP)*((s3+1)%SixP)*(s3%SixP) )%SixP )/6
                 + ( ( ((s2+1)%TwoP)*(s2%TwoP) ) % TwoP )/2 + s1%P ) % P;
    }
};

template <typename Index>
struct entry_hashfunction
{
};

template <>
struct entry_hashfunction<Index1D>
{
    entry_hashfunction(void)
    {
        for (int i=0; i<JMAX+JMINOFFSET+2; ++i) {
            power2i[i] = 1L << i;
        }
    }
    size_t power2i[JMAX+JMINOFFSET+2];

    inline
    size_t operator()(const Entry<Index1D>& entry) const
    {
        size_t l1 = power2i[entry.col_index.j+entry.col_index.xtype+JMINOFFSET] + entry.col_index.k;
        size_t l2 = power2i[entry.row_index.j+entry.row_index.xtype+JMINOFFSET] + entry.row_index.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t P=SIZELARGEHASHINDEX2D, twoP=2*SIZELARGEHASHINDEX2D;
        return (((((s2+1)%(twoP)) * (s2 % twoP)) % twoP)/2 + s1 % P) % P;
    }
    /*
    inline
    size_t operator()(const Entry<Index1D>& entry) const
    {
        int val1 = entry.col_index.xtype;
        int val2 = entry.row_index.xtype;
        val1 = (((val1 << 16) | (unsigned short) entry.col_index.j));
        val2 = (((val2 << 16) | (unsigned short) entry.row_index.j));

        std::size_t hash_value = 0;
        boost::hash_combine(hash_value, val1);
        boost::hash_combine(hash_value, val2);
        boost::hash_combine(hash_value, entry.col_index.k);
        boost::hash_combine(hash_value, entry.row_index.k);

        return hash_value;
    }
    */
};

} //namespace lawa

#include <lawa/methods/adaptive/datastructures/index.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H
