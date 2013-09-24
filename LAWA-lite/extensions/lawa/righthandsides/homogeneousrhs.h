#ifndef LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H
#define LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H 1

namespace lawa{

template <typename T>
struct HomogeneousRHS{
    
    T 
    operator()(T time, XType xtype, int j, long k) const;
    
};

} // namespace lawa

#include <lawa/righthandsides/homogeneousrhs.tcc>

#endif // LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H

