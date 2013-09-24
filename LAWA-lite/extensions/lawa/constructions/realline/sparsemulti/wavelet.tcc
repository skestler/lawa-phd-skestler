#include <cassert>
#include <iostream>
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_wavelet_evaluator.h>

namespace lawa {
  

template <typename T>
Wavelet<T,Primal,R,SparseMulti>::Wavelet(int _d)
    : d(_d), vanishingMoments(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 4: _numSplines = 4;

                _evaluator = new Evaluator[4];
                _evaluator[0] = _sparsemulti_cubic_wavelet_inner_evaluator0;
                _evaluator[1] = _sparsemulti_cubic_wavelet_inner_evaluator1;
                _evaluator[2] = _sparsemulti_cubic_wavelet_inner_evaluator2;
                _evaluator[3] = _sparsemulti_cubic_wavelet_inner_evaluator3;
            
                _support = new Support<T>[4];
                _support[0] = Support<T>(  0.,2.);
                _support[1] = Support<T>(  0.,2.);
                _support[2] = Support<T>( -2.,2.);
                _support[3] = Support<T>( -2.,2.);
            
                _singularSupport = new DenseVector<Array<T> >[4];
                _singularSupport[0] = linspace(0.0,2.0,5);
                _singularSupport[1] = linspace(0.0,2.0,5);
                _singularSupport[2] = linspace(-2.0,2.0,9);
                _singularSupport[3] = linspace(-2.0,2.0,9);

                _ScalingFactors.engine().resize(4,0);
                _ScalingFactors = 15./std::sqrt(13./7.), 39./(4.*std::sqrt(11./7.)),
                                  60./std::sqrt(2467613./2002.), 60./std::sqrt(7841./2002.);

                _max_support = Support<T>(-2.,2.);
                break;

        default: std::cerr << "Wavelet<T,Primal,R,SparseMulti> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}

template <typename T>
Wavelet<T,Primal,R,SparseMulti>::Wavelet(const Basis<T,Primal,R,SparseMulti> &_basis)
    : d(_basis.d), vanishingMoments(_basis.d)
{
    assert(d>=2);

    switch (d) {
        case 4: _numSplines = 4;

                _evaluator = new Evaluator[4];
                _evaluator[0] = _sparsemulti_cubic_wavelet_inner_evaluator0;
                _evaluator[1] = _sparsemulti_cubic_wavelet_inner_evaluator1;
                _evaluator[2] = _sparsemulti_cubic_wavelet_inner_evaluator2;
                _evaluator[3] = _sparsemulti_cubic_wavelet_inner_evaluator3;

                _support = new Support<T>[4];
                _support[0] = Support<T>(  0.,2.);
                _support[1] = Support<T>(  0.,2.);
                _support[2] = Support<T>( -2.,2.);
                _support[3] = Support<T>( -2.,2.);

                _singularSupport = new DenseVector<Array<T> >[4];
                _singularSupport[0] = linspace(0.0,2.0,5);
                _singularSupport[1] = linspace(0.0,2.0,5);
                _singularSupport[2] = linspace(-2.0,2.0,9);
                _singularSupport[3] = linspace(-2.0,2.0,9);

                _ScalingFactors.engine().resize(4,0);
                _ScalingFactors = 15./(2.*std::sqrt(13./7.)), 39./(8.*std::sqrt(11./7.)),
                                  30./std::sqrt(2467613./2002.), 30./std::sqrt(7841./2002.);
//                _ScalingFactors = 15./std::sqrt(13./7.), 39./(4.*std::sqrt(11./7.)),
//                                  60./std::sqrt(2467613./2002.), 60./std::sqrt(7841./2002.);

                _max_support = Support<T>(-2.,2.);
                break;

        default: std::cerr << "Wavelet<T,Primal,R,SparseMulti> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}
    
template <typename T>
Wavelet<T,Primal,R,SparseMulti>::~Wavelet()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
Wavelet<T,Primal,R,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);

    return pow2ih<T>(2*j*deriv+j) * _ScalingFactors(type) *
    _evaluator[type](pow2i<T>(j)*x - shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Primal,R,SparseMulti>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
//    std::cerr << "k = " << k << " : type = " << type << ", shift = " << shift << std::endl;

    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
Support<T>
Wavelet<T,Primal,R,SparseMulti>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,SparseMulti>::singularSupport(int j, long k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;
}

template <typename T>
T
Wavelet<T,Primal,R,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

template <typename T>
long
Wavelet<T,Primal,R,SparseMulti>::_shift(long k) const
{
    if (d==4) {
        return k>=0 ? 2*( k/_numSplines ) : 2*( -((-k-1)/_numSplines+1) );
    }
    else {
        std::cerr << "Wavelet<T,Primal,R,SparseMulti> not implemented for d=" << d << std::endl;
        exit(1);
        return 1;
    }
}

template <typename T>
int
Wavelet<T,Primal,R,SparseMulti>::_type(long k) const
{
    if (d==4) {
        return k>=0 ? (int) (k%_numSplines) : (int) _numSplines-(int)((-k+_numSplines-1)%_numSplines)-1;
    }
    else {
        std::cerr << "Wavelet<T,Primal,R,SparseMulti> not implemented for d=" << d << std::endl;
        exit(1);
        return 1;
    }
}

} // namespace lawa
