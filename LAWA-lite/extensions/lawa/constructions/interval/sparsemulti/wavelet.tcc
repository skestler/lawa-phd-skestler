#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Primal,Interval,SparseMulti>::Wavelet(const Basis<T,Primal,Interval,SparseMulti> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
}
    
template <typename T>
Wavelet<T,Primal,Interval,SparseMulti>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Primal,Interval,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2ih<T>(2*j*deriv+j) * basis._leftScalingFactors(k) *
                   basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
        }
        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            int type  = (int)(k % basis._numInnerParts);
            long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            //std::cerr << "type = " << type << ", shift = " << shift << " "  << std::endl;
            return pow2ih<T>(2*j*deriv+j) * basis._innerScalingFactors(type) *
                   basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
        }

        // right boundary
        long shift = basis.cardJ(j)/2;
        return pow2ih<T>(2*j*deriv+j) * basis._rightScalingFactors(0) *
               basis._rightEvaluator[0](pow2i<T>(j)*x-shift,deriv);
    }
}
    
template <typename T>
Support<T>
Wavelet<T,Primal,Interval,SparseMulti>::support(int j, long k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSupport[k];
        }

        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            int type  = (int)(k % basis._numInnerParts);
            long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
        }

        // right boundary
        long shift = basis.cardJ(j)/2;
        return pow2i<T>(-j) * (basis._rightSupport[0]+shift);
    }
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,Interval,SparseMulti>::singularSupport(int j, long k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSingularSupport[k];
        }

        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            int type  = (int)(k % basis._numInnerParts);
            long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            DenseVector<Array<T> > result = basis._innerSingularSupport[type];
            result += shift;
            return pow2i<T>(-j) * result;
        }

        // right boundary
        long shift = basis.cardJ(j)/2;
        DenseVector<Array<T> > result = basis._rightSingularSupport[0];
        result += shift;
        return pow2i<T>(-j) * result;
    }
}
    
template <typename T>
T
Wavelet<T,Primal,Interval,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}
    
} // namespace lawa
