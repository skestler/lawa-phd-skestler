#ifndef LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H
#define LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H 1

#include <cassert>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.h>

namespace lawa {
    
template <typename _T>
class MRA<_T,Primal,RPlus,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = RPlus;
        static const Construction Cons = SparseMulti;
        
        typedef BasisFunction<T,Primal,RPlus,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,RPlus,SparseMulti>       BSplineType;
        
        MRA(int d, int j=1);
        
        ~MRA();
        
        // cardinalities of left index sets.
                long
        cardIL(int j=1) const;
        
        // ranges of left index sets.
        Range<long>
        rangeIL(int j=-1) const;
        
        int
        level() const;
        
        void
        setLevel(int j) const;
        
        template <BoundaryCondition BC>
        void
        enforceBoundaryCondition();
        
        const int d;     
        const int j0;          // minimal used(!) level.
        
        BSpline<T,Primal,RPlus,SparseMulti> phi;
        
    private:
        DenseVector<Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
        
        mutable int _j;                // the current level.
    
        friend class BSpline<T,Primal,RPlus,SparseMulti>;
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        unsigned int _numLeftParts,
                     _numInnerParts;
        Evaluator *_leftEvaluator,
                  *_innerEvaluator;
        Support<T> *_leftSupport,
                   *_innerSupport;
        DenseVector<Array<T> > *_leftSingularSupport,
                               *_innerSingularSupport;
        DenseVector<Array<T> > _leftScalingFactors,
                               _innerScalingFactors;
};
    
} // namespace lawa

#include <lawa/constructions/rplus/sparsemulti/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H
