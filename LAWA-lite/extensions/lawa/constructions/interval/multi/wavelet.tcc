#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC 1

#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::Wavelet(const Basis<T,Orthogonal,Interval,Multi> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
    switch (d) {
        case 1:
            initialticsize = pow2i<T>(-1);
            break;
        case 2:
            initialticsize = pow2i<T>(-3);
            break;

        case 3:
            initialticsize = pow2i<T>(-4);
            break;

        case 4:
            initialticsize = pow2i<T>(-4);
            break;

        default: std::cerr << "Wavelet<T,Orthogonal,Interval,Multi> not yet realized"
                    " for d = " << d << ". Stopping..." << std::endl;
                    exit(1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    k -= 1;
    // left boundary
    if (k<basis._numLeftParts) {
        //std::cerr << " Left boundary: k = " << k << std::endl;
        return pow2ih<T>(2*j*deriv+j) * basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //std::cerr << "  k = " << k << " : type = " << type << ", shift = " << shift << std::endl;
        return pow2ih<T>(2*j*deriv+j) * 
        basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * basis._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,Interval,Multi>::support(int j, long k) const
{
    k -= 1;
    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j) -1 - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
    return pow2i<T>(-j) * (basis._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,Interval,Multi>::singularSupport(int j, long k) const
{
    k -= 1;    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSingularSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        DenseVector<Array<T> > result = basis._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j)-1 - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1. - basis._numLeftParts)/basis._numInnerParts);
    DenseVector<Array<T> > result = basis._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}
    
template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    //todo: Critical! Yields totally wrong results if too small, increases run-time when too large!
    //return pow2i<T>(-(j+4));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
Wavelet<T,Orthogonal,Interval,Multi>::getRefinement(int j, long k,
                                                    int &refinement_j, long &refinement_k_first) const
{
    k -= 1;
    refinement_j = j + basis._addRefinementLevel;
    // left boundary
    if (k<basis._numLeftParts) {
        refinement_k_first = basis._leftOffsets[k];
        return &(basis._leftRefCoeffs[k]);
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = (long)lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //refinement_k_first = pow2i<long>(basis._addRefinementLevel)*shift+basis._innerOffsets[type];
        refinement_k_first = basis._shiftFactor*shift+basis._innerOffsets[type];
        return &(basis._innerRefCoeffs[type]);
    }
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    long shift = (long)pow2i<long>(j)-1;
    //refinement_k_first = pow2i<long>(basis._addRefinementLevel)*shift+basis._rightOffsets[type];
    refinement_k_first = basis._shiftFactor*shift+basis._rightOffsets[type];
    return &(basis._rightRefCoeffs[type]);
}

template <typename T>
int
Wavelet<T,Orthogonal,Interval,Multi>::getRefinementLevel(int j) const
{
    return j + basis._addRefinementLevel;
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::getL2Norm(int j, long k) const
{
    return 1.;
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::getH1SemiNorm(int j, long k) const
{
    long double pow2ij = (long double)(1L << j);
    k -= 1;
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2ij * basis._leftH1SemiNorms[k];
    }

    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        return pow2ij * basis._innerH1SemiNorms[type];
    }

    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    return pow2ij * basis._rightH1SemiNorms[type];
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
