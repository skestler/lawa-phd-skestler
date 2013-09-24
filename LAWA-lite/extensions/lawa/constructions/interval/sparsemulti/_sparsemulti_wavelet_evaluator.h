#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_WAVELET_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_WAVELET_EVALUATOR_H 1

#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.h>

namespace lawa {

//--- cubic evaluators --------------------------------------------------------
    
template <typename T>
    T
    _sparsemulti_cubic_wavelet_inner_evaluator0(T x, unsigned short deriv);
    
template <typename T>
    T
    _sparsemulti_cubic_wavelet_inner_evaluator1(T x, unsigned short deriv);
    
template <typename T>
    T
    _sparsemulti_cubic_wavelet_inner_evaluator2(T x, unsigned short deriv);

template <typename T>
    T
    _sparsemulti_cubic_wavelet_inner_evaluator3(T x, unsigned short deriv);

template <typename T>
    T
    _sparsemulti_cubic_wavelet_left_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _sparsemulti_cubic_wavelet_left_evaluator1(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/_sparsemulti_wavelet_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_WAVELET_EVALUATOR_H
