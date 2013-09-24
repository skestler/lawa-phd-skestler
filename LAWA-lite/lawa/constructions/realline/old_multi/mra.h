#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_MRA_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_MRA_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>

namespace lawa {

template <typename _T>
class MRA<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;

        typedef BasisFunction<T,Orthogonal,R,Multi> BasisFunctionType;
        typedef BSpline<T,Orthogonal,R,Multi>       BSplineType;

        MRA(int _d, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        const int d;
        const int j0;          // minimal used(!) level.

        BSpline<T,Orthogonal,R,Multi> phi;

    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/realline/multi/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_REALLINE_MRA_H

