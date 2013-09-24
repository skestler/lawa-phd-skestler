#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
BSpline<T,Primal,R,SparseMulti>::BSpline(const int _d)
    : d(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 4: _numSplines = 2;
                _evaluator = new Evaluator[2];
                _evaluator[0] = _cubic_sparsemulti_scaling_inner_evaluator0;
                _evaluator[1] = _cubic_sparsemulti_scaling_inner_evaluator1;
            
                _support = new Support<T>[2];
                _support[0] = Support<T>(-1.,1.);
                _support[1] = Support<T>(-1.,1.);

                _singularSupport = new DenseVector<Array<T> >[2];
                _singularSupport[0] = linspace(-1.,1.,3);
                _singularSupport[1] = linspace(-1.,1.,3);

                _ScalingFactors.engine().resize(2,0);
                _ScalingFactors = std::sqrt(35./26.),std::sqrt(105./2.);

                _max_support = Support<T>(-1.,1.);
                break;
            
        default: std::cerr << "BSpline<T,Primal,R,SparseMulti> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }

}
    
template <typename T>
BSpline<T,Primal,R,SparseMulti>::~BSpline()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
BSpline<T,Primal,R,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);

    return pow2ih<T>(2*j*deriv+j) * _ScalingFactors(type) *
           _evaluator[type](pow2i<T>(j)*x - shift, deriv);
                                        
}
    
template <typename T>
Support<T>
BSpline<T,Primal,R,SparseMulti>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
Support<T>
BSpline<T,Primal,R,SparseMulti>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Primal,R,SparseMulti>::singularSupport(int j, long k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;    
}
    
template <typename T>
T
BSpline<T,Primal,R,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

template <typename T>
long
BSpline<T,Primal,R,SparseMulti>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
BSpline<T,Primal,R,SparseMulti>::_type(long k) const
{
    if (d==4) {
        return k>=0 ? (int) k%_numSplines : (int) _numSplines - (int)(-k+1)%_numSplines - 1;
    }
    else {
        std::cerr << "BSpline<T,Primal,R,SparseMulti> not implemented for d=" << d << std::endl;
        exit(1);
        return 1;
    }
}

    
} // namespace lawa
