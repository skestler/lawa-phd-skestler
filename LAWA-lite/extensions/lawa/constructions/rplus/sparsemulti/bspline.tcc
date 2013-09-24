#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
BSpline<T,Primal,RPlus,SparseMulti>::BSpline(const MRA<T,Primal,RPlus,SparseMulti> &_mra)
    : mra(_mra), d(_mra.d)
{
    switch (d) {
            case 4: _numSplines = 2;
                    _max_support = Support<T>(0,2);
                    break;

            default: std::cerr << "BSpline<T,Primal,RPlus,SparseMulti> not yet realized"
                                          " for d = " << d << ". Stopping." << std::endl;
                     exit(-1);
    }
}
    
template <typename T>
BSpline<T,Primal,RPlus,SparseMulti>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,RPlus,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(0<=k);
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2ih<T>(2*(j)*deriv+j) * mra._leftScalingFactors(0) *
                   mra._leftEvaluator[0](pow2i<T>(j)*x, deriv);
        }

         // inner part
        int type  = (int)(k % mra._numInnerParts);
        long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
        //std::cerr << "type = " << type << ", shift = " << shift << " " << mra._innerScalingFactors(type) << std::endl;
        return pow2ih<T>(2*j*deriv+j) * mra._innerScalingFactors(type) *
               mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
}
    
template <typename T>
Support<T>
BSpline<T,Primal,RPlus,SparseMulti>::support(int j, long k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSupport[0];
        }
        // inner part
        int type  = (int)(k % mra._numInnerParts);
        long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
        return pow2i<T>(-j) * (mra._innerSupport[type]+shift);

    }
}

template <typename T>
Support<T>
BSpline<T,Primal,RPlus,SparseMulti>::max_support() const
{
    return _max_support;
}


template <typename T>
DenseVector<Array<T> >
BSpline<T,Primal,RPlus,SparseMulti>::singularSupport(int j, long k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSingularSupport[0];
        }

        int type  = (int)(k % mra._numInnerParts);
        long shift = (long)std::ceil(T(k) / T(mra._numInnerParts));
        DenseVector<Array<T> > result = mra._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
}


template <typename T>
T
BSpline<T,Primal,RPlus,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

        
} // namespace lawa
