#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H 1

namespace lawa {

//--- cubic evaluators -----------------------------------------------------

template <typename T>
    T
    _cubic_sparsemulti_scaling_left_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_right_evaluator0(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H
