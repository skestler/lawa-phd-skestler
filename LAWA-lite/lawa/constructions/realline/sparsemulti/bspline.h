#ifndef LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Primal,R,SparseMulti>
    : public BasisFunction<_T,Primal,R,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = R;
        static const Construction Cons = SparseMulti;
        
        BSpline(const int _d);

        //TODO    BSpline(MRA<T,Primal,R,SparseMulti> &mra);
        
        virtual
        ~BSpline();
        
        T
        operator()(T x, int j, long k, unsigned short deriv) const;
        
        Support<T>
        support(int j, long k) const;
        
        Support<T>
        max_support() const;

        DenseVector<Array<T> >
        singularSupport(int j, long k) const;

        T
        tic(int j) const;

        const unsigned int d;
        unsigned int _numSplines;

    private:
        typedef T (*Evaluator)(T x, unsigned short deriv);

        long
        _shift(long k) const;

        int
        _type(long k) const;

        Evaluator *_evaluator;
        Support<T> *_support;
        DenseVector<Array<T> > *_singularSupport;
        DenseVector<Array<T> >  _ScalingFactors;

        Support<T> _max_support;
    //        T
    //TODO    tic(int j) const;
    //    int polynomialOrder;
};

} // namespace lawa

#include <lawa/constructions/realline/sparsemulti/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_BSPLINE_H
