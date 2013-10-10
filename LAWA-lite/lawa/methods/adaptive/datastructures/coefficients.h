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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H 1

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_map.hpp>
#else
    #include <ext/hash_map>
#endif

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

template <SortingCriterion S, typename T, typename Index>
struct Coefficients
{
};

/* **********************************************************************************************
 * Lexicographically sorted coefficient vector
 * ********************************************************************************************** */

template <typename T, typename Index>
#ifdef TRONE
struct Coefficients<Lexicographical,T,Index> : public std::tr1::unordered_map<Index, T,
                                                                       index_hashfunction<Index>,
                                                                       index_eqfunction<Index> >
#elif BOOST
struct Coefficients<Lexicographical,T,Index> : public boost::unordered_map<Index, T,
                                                                       index_hashfunction<Index>,
                                                                       index_eqfunction<Index> >
#else
struct Coefficients<Lexicographical,T,Index> : public __gnu_cxx::hash_map<Index, T,
                                                                       index_hashfunction<Index>,
                                                                       index_eqfunction<Index> >
#endif
{

    Coefficients(void); //required in rhs.h

    #ifdef TRONE
    Coefficients(size_t n)
     :std::tr1::unordered_map<Index, T, index_hashfunction<Index>, index_eqfunction<Index> >::unordered_map(n)
     {

     }
    void
    Rehash(size_t n) { this->rehash(n); }
    #elif BOOST
    Coefficients(size_t n)
     :boost::unordered_map<Index, T, index_hashfunction<Index>, index_eqfunction<Index> >::unordered_map(n)
     {

     }
    void
    Rehash(size_t n) { this->rehash(n); }
    #else
    Coefficients(size_t n)
    :__gnu_cxx::hash_map<Index, T, index_hashfunction<Index>, index_eqfunction<Index> >::hash_map(n) {

    };
    void
    Rehash(size_t n) { this->resize(n); }
    #endif




    Coefficients<Lexicographical,T,Index>&
    operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index>&
    operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index>
    operator-(const Coefficients<Lexicographical,T,Index> &_coeff) const;

    Coefficients<Lexicographical,T,Index> &
    operator-=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index> &
    operator+=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index> &
    operator*=(const T factor);

    Coefficients<Lexicographical,T,Index>
    operator+(const Coefficients<Lexicographical,T,Index> &_coeff) const;

    T
    operator*(const Coefficients<Lexicographical,T,Index> &_coeff) const;
    
    //todo:: revise!!!
    Coefficients<Lexicographical,T,Index>
    operator*(const T factor) const;
    
    void
    scale(const T factor);

    void
    setToZero();

    T
    norm(T tau=2.0) const;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Lexicographical,T,Index> &c);

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
operator*(T alpha, const Coefficients<Lexicographical,T,Index> &_coeff);

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
P(const Coefficients<Lexicographical,T,Index> &v, const IndexSet<Index> &Lambda);

template <typename T, typename Index>
void
P(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &v);

template <typename T, typename Index>
IndexSet<Index>
supp(const Coefficients<Lexicographical,T,Index> &v);

template <typename T, typename Index>
void
FillWithZeros(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &_coeff);

template <typename T>
void
getLevelInfo(const Coefficients<Lexicographical,T,Index2D> &Lambda, Index2D &maxIndex,
             Index2D &maxWaveletIndex, int *jmax, int &arrayLength);

template <typename T>
void
getLevelInfo(const Coefficients<Lexicographical,T,Index3D> &Lambda, Index3D &maxIndex,
             Index3D &maxWaveletIndex, int *jmax, int &arrayLength);

template <typename T, typename Index>
void
Pe1(const Coefficients<Lexicographical,T,Index> &v,
    Coefficients<Lexicographical,T,Index1D> &Pe1_v);

template <typename T, typename Index>
void
Pe2(const Coefficients<Lexicographical,T,Index> &v,
    Coefficients<Lexicographical,T,Index1D> &Pe2_v);

template <typename T, typename Index>
void
Pe3(const Coefficients<Lexicographical,T,Index> &v,
    Coefficients<Lexicographical,T,Index1D> &Pe3_v);


/* **********************************************************************************************
 * Bucket sorted coefficient vector
 * ********************************************************************************************** */


template <typename T, typename Index>
struct Coefficients<Bucket,T,Index>
{

    typedef typename std::list<const std::pair<const Index,T>* >               BucketEntry;
    typedef typename std::vector<BucketEntry> Buckets;


    Coefficients();        //required in rhs.h

    void
    bucketsort(const Coefficients<Lexicographical,T,Index> &_coeff, T eps);

    int
    addBucketToIndexSet(IndexSet<Index> &Lambda, int bucketnumber);

    void
    addBucketToCoefficients(Coefficients<Lexicographical,T,Index> &coeff, int bucketnumber);


    T supremumnorm;

    //std::vector<Coefficients<Lexicographical,T,Index> > buckets;
    Buckets                    buckets;
    std::vector<long double>   bucket_ell2norms;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Bucket,T,Index> &c);


/* **********************************************************************************************
 * Coefficient vector sorted by absolute values
 * ********************************************************************************************** */

template <typename T, typename Index>
struct Coefficients<AbsoluteValue,T,Index> : std::multimap<T,Index,lt<AbsoluteValue,T> >
{
    using std::multimap<T,Index,lt<AbsoluteValue,T> >::insert;
    using std::multimap<T,Index,lt<AbsoluteValue,T> >::erase;
    
    Coefficients();

    Coefficients<AbsoluteValue,T,Index>&
    operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<AbsoluteValue,T,Index>&
    operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

    T
    norm(T tau=2.0) const;

    T
    l2bestnterm(int n) const;

    T
    wtauNorm(T tau) const;

    DenseVector<Array<T> >
    norm_sections() const;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<AbsoluteValue,T,Index> &c);

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/coefficients.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H

//struct Coefficients<Lexicographical,T,Index> : std::map<Index,T,lt<Lexicographical,Index> >
//using std::map<Index,T,lt<Lexicographical,Index> >::insert;
//using std::map<Index,T,lt<Lexicographical,Index> >::erase;
//using __gnu_cxx::hash_map<Index, T, index_hashfunction<Index>, index_eqfunction<Index> >::hash_map;

