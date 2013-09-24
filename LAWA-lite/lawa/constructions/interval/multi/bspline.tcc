#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_TCC 1

#include <cassert>
#include <iostream>
#include <lawa/math/math.h>

namespace lawa {

template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::BSpline(const MRA<T,Orthogonal,Interval,Multi> &_mra)
    : mra(_mra), d(_mra.d), initialticsize(pow2i<T>(-3))
{
    switch (d) {
        case 1:
            initialticsize = 1.;
            break;
        case 2:
            initialticsize = pow2i<T>(-2);
            break;

        case 3:
            initialticsize = pow2i<T>(-3);
            break;

        case 4:
            initialticsize = pow2i<T>(-2);
            break;

        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
                    " for d = " << d << ". Stopping." << std::endl;
                    exit(-1);
    }

}
    
template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2ih<T>(2*j*deriv+j) * mra._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2ih<T>(2*j*deriv+j) * 
               mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    // right part
    int type  = (int)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    return pow2ih<T>(2*j*deriv+j) * mra._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
BSpline<T,Orthogonal,Interval,Multi>::support(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSupport[k];
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
    }
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    return pow2i<T>(-j) * (mra._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Orthogonal,Interval,Multi>::singularSupport(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSingularSupport[k];
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        DenseVector<Array<T> > result = mra._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    DenseVector<Array<T> > result = mra._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    //return pow2i<T>(-(j+3));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
BSpline<T,Orthogonal,Interval,Multi>::getRefinement(int j, long k,
                                                    int &refinement_j, long &refinement_k_first) const
{
    refinement_j = j + mra._addRefinementLevel;
    // left boundary
    if (k<mra._numLeftParts) {
        refinement_k_first = mra._leftOffsets[k];
        return &(mra._leftRefCoeffs[k]);
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = (long)iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        refinement_k_first = mra._shiftFactor*shift+mra._innerOffsets[type];
        return &(mra._innerRefCoeffs[type]);
    }
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    refinement_k_first = mra._shiftFactor*shift+mra._rightOffsets[type];
    return &(mra._rightRefCoeffs[type]);
}

template <typename T>
int
BSpline<T,Orthogonal,Interval,Multi>::getRefinementLevel(int j) const
{
    return j + mra._addRefinementLevel;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::getL2Norm(int j, long k) const
{
    return 1.;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::getH1SemiNorm(int j, long k) const
{
    long double pow2ij = (long double)(1L << j);
    if (k<mra._numLeftParts) {
        return pow2ij * mra._leftH1SemiNorms[k];
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        return pow2ij * mra._innerH1SemiNorms[type];
    }
    // right part
    int type  = (int)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    return pow2ij * mra._rightH1SemiNorms[type];
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BSPLINE_TCC
