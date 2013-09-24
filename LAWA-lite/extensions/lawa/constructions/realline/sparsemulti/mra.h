#ifndef LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H
#define LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>

namespace lawa {

template <typename _T>
class MRA<_T,Primal,R,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = R;
        static const Construction Cons = SparseMulti;

        typedef BasisFunction<T,Primal,R,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,R,SparseMulti>       BSplineType;

        MRA(int _d, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        const int d;
        const int j0;          // minimal used(!) level.

        BSpline<T,Primal,R,SparseMulti> phi;

    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/realline/sparsemulti/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H

