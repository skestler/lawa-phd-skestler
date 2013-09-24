#ifndef LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_H
#define LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_H 1

#include <cassert>
#include <list>
#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T>
    T
    w(int i, int d, const DenseVector<Array<T> > &knots, T x);

template <typename T>
    GeMatrix<FullStorage<T,cxxblas::ColMajor> >
    insertKnot(int d, DenseVector<Array<T> > &knots, T x);

} // namespace lawa

#include <lawa/constructions/interval/primbs/primal/spline_helper.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_H

