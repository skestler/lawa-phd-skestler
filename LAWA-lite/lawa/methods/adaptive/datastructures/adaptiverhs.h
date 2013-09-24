#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ADAPTIVERHS_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ADAPTIVERHS_H 1

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {
    
template <typename T, typename Index>
struct AdaptiveRhs {
               
    virtual T
    operator()(const Index &lambda) = 0;

    virtual Coefficients<Lexicographical,T,Index>
    operator()(const IndexSet<Index> &Lambda) = 0;

    virtual Coefficients<Lexicographical,T,Index>
    operator()(T tol) = 0;

        
};
    
} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_ADAPTIVERHS_H
