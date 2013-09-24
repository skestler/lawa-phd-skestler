#ifndef LAWA_OPERATORS_UNIFORMOPERATOR2D_H
#define LAWA_OPERATORS_UNIFORMOPERATOR2D_H 1

#include <lawa/operators/operator2d.h>

namespace lawa {
    
template <typename T>
struct UniformOperator2D : public Operator2D<T> {
	
    virtual T
    operator()(XType row_xtype_x, int j1_x, long k1_x,
               XType row_xtype_y, int j1_y, long k1_y,
               XType col_xtype_x, int j2_x, long k2_x,
               XType col_xtpye_y, int j2_y, long k2_y) const = 0;
    
    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) const = 0;
    
};
    
} // namespace lawa

#endif // LAWA_OPERATORS_UNIFORMOPERATOR2D_H
