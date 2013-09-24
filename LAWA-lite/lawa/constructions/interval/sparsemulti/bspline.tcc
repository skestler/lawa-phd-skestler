#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
BSpline<T,Primal,Interval,SparseMulti>::BSpline(const MRA<T,Primal,Interval,SparseMulti> &_mra)
    : mra(_mra), d(_mra.d)
{
}
    
template <typename T>
BSpline<T,Primal,Interval,SparseMulti>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,Interval,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(0<=k && k<mra.cardI(j));
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2ih<T>(2*(j)*deriv+j) * mra._leftScalingFactors(0) *
                   mra._leftEvaluator[0](pow2i<T>(j)*x, deriv);
        }
        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            int type  = (int)(k % mra._numInnerParts);
            long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
            //std::cerr << "type = " << type << ", shift = " << shift << " " << mra._innerScalingFactors(type) << std::endl;
            return pow2ih<T>(2*j*deriv+j) * mra._innerScalingFactors(type) *
                   mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
        }
        // right boundary
        long shift = k / mra._numInnerParts + 1;
        return pow2ih<T>(2*j*deriv+j) * mra._rightScalingFactors(0) *
               mra._rightEvaluator[0](pow2i<T>(j)*x-shift,deriv);
    }
}
    
template <typename T>
Support<T>
BSpline<T,Primal,Interval,SparseMulti>::support(int j, long k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSupport[0];
        }

        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            int type  = (int)(k % mra._numInnerParts);
            long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
            return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
        }

        // right boundary
        long shift = k / mra._numInnerParts + 1;
        return pow2i<T>(-j) * (mra._rightSupport[0]+shift);
    }
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Primal,Interval,SparseMulti>::singularSupport(int j, long k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSingularSupport[0];
        }

        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            int type  = (int)(k % mra._numInnerParts);
            long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
            DenseVector<Array<T> > result = mra._innerSingularSupport[type];
            result += shift;
            return pow2i<T>(-j) * result;
        }

        // right part
        long shift = k / mra._numInnerParts + 1;
        DenseVector<Array<T> > result = mra._rightSingularSupport[0];
        result += shift;
        return pow2i<T>(-j) * result;
    }
}

template <typename T>
T
BSpline<T,Primal,Interval,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}
        
} // namespace lawa
