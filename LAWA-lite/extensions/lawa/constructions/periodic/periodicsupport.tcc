#include <cassert>

namespace lawa{

template <typename T>
PeriodicSupport<T>::PeriodicSupport()
    : Support<T>(), li1(0), li2(0)
{
}

template <typename T>
PeriodicSupport<T>::PeriodicSupport(const T &a, const T &b)
   : Support<T>(a,b), li1(0), li2(0)
{
    assert(0 <= Support<T>::l1);
    assert(l2 <= 1.);
}

template <typename T>
PeriodicSupport<T>::PeriodicSupport(const T&a, const T &b, const T &ai, const T &bi)
   : Support<T>(a,b), li1(ai), li2(bi)
{
    assert(l1 <= l2);
    assert(0 <= l1);
    assert(l2 <= 1.);

    if(li1 != li2){
        assert(li1 < li2);
        assert(li1 >= l1);
        assert(li2 <= l2);
        assert(l1 == 0.);
        assert(l2 == 1.);
    }
    else{
        li1 = 0;
        li2 = 0;
    }
}

template <typename T>
PeriodicSupport<T>::~PeriodicSupport()
{
}

template <typename T>
T
PeriodicSupport<T>::length() const
{
    if(gaplength() > 0){
        return (li1 - l1) + (l2 -li2);
    }
    else{
        return l2-l1;
    }
}

template <typename T>
T
PeriodicSupport<T>::gaplength() const
{
    return li2 - li1;
}

//------------------------------------------------------------------------------------

template <typename T>
bool
inner(T x, const PeriodicSupport<T> &supp) {
    if(supp.gaplength() > 0){
        return (x >= supp.l1) && (x <= supp.li1) && (x <= supp.l2) && (x >= supp.li2);
    }
    else{
        return (x >= supp.l1) && (x <= supp.l2);
    }
}

template <typename T>
T
overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2)
{
    // ''Normal'' supports:
    if((supp1.gaplength() == 0) && (supp2.gaplength() == 0)){
        return std::min(supp1.l2, supp2.l2) - std::max(supp1.l1, supp2.l1);
    }
    // Two divided supports:
    if((supp1.gaplength() > 0) && (supp2.gaplength() > 0)){
        assert((supp1.l1 == 0) && (supp1.l2 == 1));
        assert((supp2.l1 == 0) && (supp2.l2 == 1));
        return 1;
    }
    // Supp1 is divided, Supp2 not:
    if(supp1.gaplength() > 0){
        // Overlaps both right and left -> return interval covering both
        if((supp2.l2 - supp1.li2 > 0) && (supp1.li1 - supp2.l1 > 0)){
            return supp2.length();
        }
        // No or only one overlapping interval
        else{
            return std::max(std::max(0., supp2.l2 - supp1.li2), std::max(0., supp1.li1 - supp2.l1));
        }
    }

    // Supp2 is divided, Supp1 not:
    if((supp1.l2 - supp2.li2 > 0) && (supp2.li1 - supp1.l1 > 0)){
        return supp1.length();
    }
    else{
        return std::max(std::max(0., supp1.l2 - supp2.li2), std::max(0., supp2.li1 - supp1.l1));
    }
}

template <typename T>
T
overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2)
{
    // ''Normal'' supports:
    if(supp1.gaplength() == 0){
        return std::min(supp1.l2, supp2.l2) - std::max(supp1.l1, supp2.l1);
    }
    else{
        // Overlaps both right and left -> return interval covering both
        if( (supp2.l2 - supp1.li2 > 0) && (supp1.li1 - supp2.l1 > 0)){
            return supp2.length();
        }
        // No or only one overlapping interval
        else{
            return std::max(std::max(0., supp2.l2 - supp1.li2), std::max(0., supp1.li1 - supp2.l1));
        }
    }
}

template <typename T>
T
overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2)
{
    return overlap(supp2, supp1);
}

template <typename T>
T
overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2, Support<T> &common)
{
    if((supp1.gaplength() == 0) && (supp2.gaplength() == 0)){
        common.l1 = std::max(supp1.l1, supp2.l1);
        common.l2 = std::min(supp1.l2, supp2.l2);
        return common.l2 - common.l1;
    }
    if((supp1.gaplength() > 0) && (supp2.gaplength() > 0)){
        common.l1 = 0.;
        common.l2 = 1.;
        return 1;
    }
    if(supp1.gaplength() > 0){
        if((supp2.l2 - supp1.li2 > 0) && (supp1.li1 - supp2.l1 > 0)){
            common.l1 = supp2.l1;
            common.l2 = supp2.l2;
            return supp2.length();
        }
        else{
            if(supp2.l2 - supp1.li2 > 0){
                common.l1 = supp1.li2;
                common.l2 = supp2.l2;
                return common.l2 - common.l1;
            }

            if(supp1.li1 - supp2.l1 > 0){
                common.l1 = supp2.l1;
                common.l2 = supp1.li1;
                return common.l2 - common.l1;
            }

            common.l1 = 0;
            common.l2 = -1;
            return 0;
        }
    }

    if((supp1.l2 - supp2.li2 > 0) && (supp2.li1 - supp1.l1 > 0)){
        common.l1 = supp1.l1;
        common.l2 = supp1.l2;
        return supp1.length();
    }
    else{
        if(supp1.l2 - supp2.li2 > 0){
            common.l1 = supp2.li2;
            common.l2 = supp1.l2;
            return common.l2 - common.l1;
        }
        if(supp2.li1 - supp1.l1 > 0){
            common.l1 = supp1.l1;
            common.l2 = supp2.li1;
            return common.l2 - common.l1;
        }
        common.l1 = 0;
        common.l2 = -1;
        return 0;
    }
}

template <typename T>
T
overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2, Support<T> &common)
{
    // ''Normal'' supports:
    if(supp1.gaplength() == 0){
        common.l1 = std::max(supp1.l1, supp2.l1);
        common.l2 = std::min(supp1.l2, supp2.l2);
        return common.l2 - common.l1;
    }
    else{
        // Overlaps both right and left -> return interval covering both
        if((supp2.l2 - supp1.li2 > 0) && (supp1.li1 - supp2.l1 > 0)){
            common.l1 = supp2.l1;
            common.l2 = supp2.l2;
            return supp2.length();
        }
        // No or only one overlapping interval
        else{
            if(supp2.l2 - supp1.li2 > 0) {
                common.l1 = supp1.li2;
                common.l2 = supp2.l2;
                return common.l2 - common.l1;
            }

            if(supp1.li1 - supp2.l1 > 0) {
                common.l1 = supp2.l1;
                common.l2 = supp1.li1;
                return common.l2 - common.l1;
            }
            common.l1 = 0;
            common.l2 = -1;
            return 0;
        }
    }
}

template <typename T>
T
overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2, Support<T> &common)
{
    return overlap(supp2, supp1, common);
}

template <typename T>
T
minimal_overlap(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2)
{
    // ''Normal'' supports:
    if((supp1.gaplength() == 0) && (supp2.gaplength() == 0)){
        return std::min(supp1.l2, supp2.l2) - std::max(supp1.l1, supp2.l1);
    }
    // Two divided supports:
    if((supp1.gaplength() > 0) && (supp2.gaplength() > 0)){
        return 1 - std::max(supp1.li2, supp2.li2) + std::min(supp1.li1, supp2.li1);
    }
    // Supp1 is divided, Supp2 not:
    if(supp1.gaplength() > 0){
        return std::max(0., supp2.l2 - supp1.li2) + std::max(0., supp1.li1 - supp2.l1);
    }

    // Supp2 is divided, Supp1 not:
    return std::max(0., supp1.l2 - supp2.li2) + std::max(0., supp2.li1 - supp1.l1);
}

template <typename T>
T
minimal_overlap(const PeriodicSupport<T> &supp1, const Support<T> &supp2)
{
    // ''Normal'' supports:
    if(supp1.gaplength() == 0){
        return std::min(supp1.l2, supp2.l2) - std::max(supp1.l1, supp2.l1);
    }
    // Supp1 is divided
    else{
        return std::max(0., supp2.l2 - supp1.li2) + std::max(0., supp1.li1 - supp2.l1);
    }
}

template <typename T>
T
minimal_overlap(const Support<T> &supp1, const PeriodicSupport<T> &supp2)
{
    return minimal_overlap(supp2, supp1);
}

template <typename T>
T
distance(const PeriodicSupport<T> &supp1, const PeriodicSupport<T> &supp2)
{
    // ''Normal'' supports:
    if((supp1.gaplength() == 0) && (supp2.gaplength() == 0)){
        return std::max(supp1.l1, supp2.l1) - std::min(supp1.l2, supp2.l2);
    }
    // Two divided supports:
    if((supp1.gaplength() > 0) && (supp2.gaplength() > 0)){
        return 0;
    }
    // Supp1 is divided, Supp2 not:
    if(supp1.gaplength() > 0){
        // Overlaps right or left -> distance = 0
        if( (supp2.l2 - supp1.li2 >= 0) || (supp1.li1 - supp2.l1 >= 0)){
            return 0;
        }
        else{
            return std::min(supp2.l1 - supp1.li1, supp1.li2 - supp2.l2);
        }
    }

    // Supp2 is divided, Supp1 not:
    if( (supp1.l2 - supp2.li2 >= 0) || (supp2.li1 - supp1.l1 >= 0)){
        return 0;
    }
    else{
        return std::min(supp1.l1 - supp2.li1, supp2.li2 - supp1.l2);
    }

}

template <typename T>
T
distance(const PeriodicSupport<T> &supp1, const Support<T> &supp2)
{
    // ''Normal'' supports:
    if(supp1.gaplength() == 0){
        return std::max(supp1.l1, supp2.l1) - std::min(supp1.l2, supp2.l2);
    }
    else{
        // Overlaps right or left -> distance = 0
        if( (supp2.l2 - supp1.li2 >= 0) || (supp1.li1 - supp2.l1 >= 0)){
            return 0;
        }
        else{
            return std::min(supp2.l1 - supp1.li1, supp1.li2 - supp2.l2);
        }
    }
}

template <typename T>
T
distance(const Support<T> &supp1, const PeriodicSupport<T> &supp2)
{
    return distance(supp2, supp1);
}

template <typename T, typename S>
PeriodicSupport<T>
operator+(const PeriodicSupport<T> &supp, S shift)
{
    PeriodicSupport<T> newsupp;
    if(supp.gaplength() > 0){
        newsupp.l1 = supp.li2 + shift - ifloor(supp.li2 + shift);
        newsupp.l2 = supp.li1 + shift - ifloor(supp.li1 + shift);
    }
    else{
        newsupp.l1 = supp.l1 + shift - ifloor(supp.l1 + shift);
        newsupp.l2 = supp.l2 + shift - ifloor(supp.l2 + shift);
    }
    if(newsupp.l2 < newsupp.l1){
        newsupp.li1 = newsupp.l2;
        newsupp.li2 = newsupp.l1;
        newsupp.l1 = 0.;
        newsupp.l2 = 1.;
    }
    return newsupp;
}


template <typename T>
std::ostream &
operator<<(std::ostream &out, const PeriodicSupport<T> &supp)
{
    if(supp.gaplength() > 0){
        out << "[" << supp.l1 << "," << supp.li1 << "]" << " v "
            << "[" << supp.li2 << "," << supp.l2 << "]";
    }
    else{
        out << "[" << supp.l1 << "," << supp.l2 << "]";
    }
    return out;
}


} // namespace lawa

