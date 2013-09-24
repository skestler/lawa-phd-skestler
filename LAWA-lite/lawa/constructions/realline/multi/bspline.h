#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
        
        BSpline(const int _d);

        //TODO    BSpline(MRA<T,Orthogonal,R,Multi> &mra); 
        
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

        DenseVector<Array<long double> > *
        getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const;

        int
        getRefinementLevel(int j) const;

        //T
        //getL2Norm(int j, long k) const;

        //T
        //getH1SemiNorm(int j, long k) const;

        const unsigned int d;
        unsigned int _numSplines;

//    private:      // should be private one fine day

        typedef T (*Evaluator)(T x, unsigned short deriv);

        long
        _shift(long k) const;

        int
        _type(long k) const;

        Evaluator                        *_evaluator;
        Support<T>                       *_support;
        DenseVector<Array<T> >           *_singularSupport;
        Support<T>                       _max_support;

        DenseVector<Array<long double> > *_refCoeffs;
        long                             *_offsets;
        long double                      *_H1SemiNorms;
        T                                _initialticsize;
        int                              _addRefinementLevel;    //B-splines for refinement are needed on higher levels
        int                              _shiftFactor;           //Needed since we have multiple B-spline generators for refinement.
};

} // namespace lawa

#include <lawa/constructions/realline/multi/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
