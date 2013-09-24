#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H 1

namespace lawa {

template <typename X, typename Basis>
typename X::ElementType
evaluate(const Basis &basis, const int J_x, const int J_y, const flens::DenseVector<X> &coeffs,
         const typename X::ElementType x, const typename X::ElementType y, const int deriv_x,
         const int deriv_y);

} // namespace lawa

#include <lawa/methods/uniform/postprocessing/evaluate.tcc>


#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H

