#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Primal,Interval,SparseMulti>
    : public BasisFunction<_T,Primal,Interval,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = SparseMulti;

        BSpline(const MRA<T,Primal,Interval,SparseMulti> &mra);

        virtual
        ~BSpline();

        T
        operator()(T x, int j, long k, unsigned short deriv) const;

        Support<T>
        support(int j, long k) const;

        DenseVector<Array<T> >
        singularSupport(int j, long k) const;

        T
        tic(int j) const;

        const MRA<T,Primal,Interval,SparseMulti> &mra;
        const unsigned int d;
};

} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H
