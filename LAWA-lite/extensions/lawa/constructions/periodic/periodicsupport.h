#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PERIODICSUPPORT_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PERIODICSUPPORT_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/support.h>

namespace lawa {

using namespace flens;

template <typename T>
struct PeriodicSupport : public Support<T>
{
    PeriodicSupport();

    PeriodicSupport(const T &a, const T &b);

    PeriodicSupport(const T&a, const T &b, const T &ai, const T &bi);

    ~PeriodicSupport();

    T
    length() const;
    
    T
    gaplength() const;      // length of possible inner gap in support (= length[li1, li2])

    using Support<T>::l1;   // outer support bounds (maximal support: [0,1])
    using Support<T>::l2;
    T li1, li2;             // bounds of inner gap if support is [l1, li1] and [li2, l2]
                            //  if li1 = li2 (= 0), support is [l1, l2].

};

template <typename T>
    bool
    inner(T x, const PeriodicSupport<T> &supp);

template <typename T>
    T
    overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T>
    T
    overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2);

template <typename T>
    T
    overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T>
    T
    overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2, Support<T> &common);

template <typename T>
    T
    overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2, Support<T> &common);

template <typename T>
    T
    overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2, Support<T> &common);

template <typename T>
    T
    minimal_overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T>
    T
    minimal_overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2);

template <typename T>
    T
    minimal_overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T>
    T
    distance(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T>
    T
    distance(const PeriodicSupport<T> &supp1, const Support<T> &supp2);

template <typename T>
    T
    distance(const Support<T> &supp1, const PeriodicSupport<T> &supp2);

template <typename T, typename S>
    PeriodicSupport<T>
    operator+(const PeriodicSupport<T> &supp, S shift);

template <typename S, typename T>
    PeriodicSupport<T>
    operator*(S factor, const PeriodicSupport<T> &supp);

template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const PeriodicSupport<T> &supp);

} // namespace lawa

#include <lawa/constructions/periodic/periodicsupport.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PERIODICSUPPORT_H

