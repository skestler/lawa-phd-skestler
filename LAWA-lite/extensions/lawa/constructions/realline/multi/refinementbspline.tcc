#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBSPLINE_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBSPLINE_TCC 1

#include <cassert>
#include <iostream>
#include <lawa/math/math.h>

#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
BSpline<T,Orthogonal,R,MultiRefinement>
::BSpline(const int _d)
    : d(_d), _initialticsize(pow2i<T>(-3))
{
    if (d > 4) {
        std::cerr << "BSpline<T,Orthogonal,R,MultiRefinement> not yet implemented for d = " << d << std::endl;
        exit(1);
    }

    switch (d) {
        case 1:
            _shiftFactor   = 2;
            _numSplines    = 1;
            _initialticsize = 1.;
            _evaluator = new Evaluator[1];
            _evaluator[0] = _constant_bspline_inner_evaluator0;

            _support = new Support<T>[1];
            _support[0] = Support<T>(0.,1.);

            _singularSupport = new DenseVector<Array<T> >[1];
            _singularSupport[0].engine().resize(2,0);
            _singularSupport[0] = 0., 1.;

            _refCoeffs = new DenseVector<Array<long double> >[1];
            _refCoeffs[0].engine().resize(2,0);
            _refCoeffs[0] = 1.L, 1.L;
            _refCoeffs[0] *= std::pow(2.L,-0.5L);

            _offsets = new long[1];
            _offsets[0] =  0;

            break;

        case 2:
            _shiftFactor   = 2;
            _numSplines = 1;
            _initialticsize = 1.;

            _evaluator = new Evaluator[1];
            _evaluator[0] = _linear_refinement_inner_evaluator0;

            _support = new Support<T>[1];
            _support[0] = Support<T>(0.,2.);

            _singularSupport = new DenseVector<Array<T> >[1];
            _singularSupport[0].engine().resize(3,0);
            _singularSupport[0] = 0., 1., 2.;

            _refCoeffs = new DenseVector<Array<long double> >[1];
            _refCoeffs[0].engine().resize(3,0);
            _refCoeffs[0] = 0.5L, 1.L, 0.5L;
            _refCoeffs[0] *= std::pow(2.L,-0.5L);

            _offsets = new long[1];
            _offsets[0] =  0;
            break;

        case 3:
            _shiftFactor   = 2;
            _numSplines = 1;
            _initialticsize = pow2i<T>(-3);

            _evaluator = new Evaluator[1];
            _evaluator[0] = _quadratic_refinement_inner_evaluator0;
            _support = new Support<T>[1];
            _support[0] = Support<T>(0.,3.);
            _singularSupport = new DenseVector<Array<T> >[1];
            _singularSupport[0].engine().resize(4,0);
            _singularSupport[0] = 0., 1., 2., 3.;

            _refCoeffs = new DenseVector<Array<long double> >[1];
            _refCoeffs[0].engine().resize(4,0);
            _refCoeffs[0] = 0.25L, 0.75L, 0.75L, 0.25L;
            _refCoeffs[0] *= std::pow(2.L,-0.5L);
            _offsets = new long[1];
            _offsets[0] =  0;
            break;

        case 4:
            _shiftFactor = 4;
            _numSplines = 2;
            _initialticsize = pow2i<T>(-2);

            _evaluator = new Evaluator[2];
            _evaluator[0] = _cubic_refinement_inner_evaluator0;
            _evaluator[1] = _cubic_refinement_inner_evaluator1;
            _support = new Support<T>[2];
            _support[0] = Support<T>(0.,2.);
            _support[1] = Support<T>(0.,2.);
            _singularSupport = new DenseVector<Array<T> >[2];
            _singularSupport[0].engine().resize(3,0);
            _singularSupport[0] = 0., 1., 2.;
            _singularSupport[1].engine().resize(3,0);
            _singularSupport[1] = 0., 1., 2.;

            _refCoeffs = new DenseVector<Array<long double> >[2];
            _refCoeffs[0].engine().resize(5,0);
            _refCoeffs[0] = 0.25L, 0.625L, 0.75L, 0.25L, 0.125L;
            _refCoeffs[0] *= std::pow(2.L,-0.5L);
            _refCoeffs[1].engine().resize(5,0);
            _refCoeffs[1] = 0.125L, 0.25L, 0.75L, 0.625L, 0.25L;
            _refCoeffs[1] *= std::pow(2.L,-0.5L);
            _offsets = new long[2];
            _offsets[0] =  0;
            _offsets[1] =  1;
            break;

        default: std::cerr << "BSpline<T,Orthogonal,R,MultiRefinement> not yet realized"
                    " for d = " << d << ". Stopping." << std::endl;
                    exit(-1);
    }

}

template <typename T>
BSpline<T,Orthogonal,R,MultiRefinement>::~BSpline()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
    delete[] _refCoeffs;
    delete[] _offsets;
}

template <typename T>
T
BSpline<T,Orthogonal,R,MultiRefinement>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int  type  = _type(k);
    const long shift = _shift(k);

    return pow2ih<T>(2*j*deriv+j) * _evaluator[type](pow2i<T>(j)*x - shift, deriv);
}

template <typename T>
Support<T>
BSpline<T,Orthogonal,R,MultiRefinement>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);

    //std::cerr << "BSpline<T,Orthogonal,R,MultiRefinement>: k = " << k << ", type = " << type
    //          << ", shift = " << shift << std::endl;

    return pow2i<T>(-j) * (_support[type] + shift);
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Orthogonal,R,MultiRefinement>::singularSupport(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);

    DenseVector<Array<T> > result = _singularSupport[type];
    result += shift;

    return pow2i<T>(-j) * result;
}

template <typename T>
T
BSpline<T,Orthogonal,R,MultiRefinement>::tic(int j) const
{
    //return pow2i<T>(-(j+3));
    return _initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
BSpline<T,Orthogonal,R,MultiRefinement>::
getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const
{
    refinement_j = j + 1;
    long shift = this->_shift(k);
    int  type  = this->_type(k);

    refinement_k_first = _shiftFactor*shift+_offsets[type];
    return &(_refCoeffs[type]);
}

template <typename T>
int
BSpline<T,Orthogonal,R,MultiRefinement>::getRefinementLevel(int j) const
{
    return j + 1;
}

template <typename T>
long
BSpline<T,Orthogonal,R,MultiRefinement>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
BSpline<T,Orthogonal,R,MultiRefinement>::_type(long k) const
{
    return k>=0 ? k % _numSplines : _numSplines - (_numSplines-1-k) % _numSplines - 1;

}



} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_REALLINE_REFINEMENTBSPLINE_TCC
