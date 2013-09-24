#ifndef LAWA_OPERATORS_OPERATOR2D_H
#define LAWA_OPERATORS_OPERATOR2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {
    
template <typename T>
struct Operator2D {     
               
    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) = 0;
    
};
    
} // namespace lawa

#endif // LAWA_OPERATORS_OPERATOR2D_H
