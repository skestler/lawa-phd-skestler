#include <cassert>
#include <iostream>

namespace lawa {
    
template <typename T>
    T
    _linear_bspline_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_bspline_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _linear_bspline_evaluator2(T x, unsigned short deriv);
    
//------------------------------------------------------------------------------

template <typename T>
BSpline<T,Orthogonal,R,Multi>::BSpline(const int _d)
    : d(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 2: _numSplines = 3;
                
                _evaluator = new Evaluator[3];
                _evaluator[0] = _linear_bspline_evaluator0;
                _evaluator[1] = _linear_bspline_evaluator1;                
                _evaluator[2] = _linear_bspline_evaluator2;
            
                _support = new Support<T>[3];
                _support[0] = Support<T>(-1,1);
                _support[1] = Support<T>( 0,1);
                _support[2] = Support<T>( 0,1);

                _singularSupport = new DenseVector<Array<T> >[3];
                _singularSupport[0] = linspace(-1.0,1.0,9);
                _singularSupport[1] = linspace(0.0,1.0,3);
                _singularSupport[2] = linspace(0.0,1.0,5);

                _max_support = Support<T>(-1,1);
                break;
            
        default: std::cerr << "BSpline<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }

}
    
template <typename T>
BSpline<T,Orthogonal,R,Multi>::~BSpline()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
BSpline<T,Orthogonal,R,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);

    return pow2ih<T>(2*j*deriv+j) *
           _evaluator[type](pow2i<T>(j)*x - shift, deriv);
                                        
}
    
template <typename T>
Support<T>
BSpline<T,Orthogonal,R,Multi>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
Support<T>
BSpline<T,Orthogonal,R,Multi>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Orthogonal,R,Multi>::singularSupport(int j, long k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;    
}
    
template <typename T>
T
BSpline<T,Orthogonal,R,Multi>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

template <typename T>
long
BSpline<T,Orthogonal,R,Multi>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
BSpline<T,Orthogonal,R,Multi>::_type(long k) const
{
    return k>=0 ? (int) k%3 : (int) _numSplines - (int)(-k+2)%((int)_numSplines) - 1;
}

//------------------------------------------------------------------------------
    
template <typename T>
T
_linear_bspline_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0 <= x && x < 0.25) {
            value = 2.171240593367237661667 - 11.00419101450331474149* x;
        } else if (-0.75 <= x &&  x < -0.5){
            value = -2.458152258280902698695 - 3.441413161593263778173* x;
        } else if (0.5 <= x && x < 0.75){
            value = 0.484297173401595733543 - 0.678016042762234026960* x;
        } else if (0.75 <= x && x < 1.0){
            value = -0.096859434680319146709 + 0.096859434680319146709* x;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 + 0.491630451656180539739* x;
        } else if ( 0.25 <= x && x < 0.5){
            value = -1.30490347253766076747 + 2.90038524911627897507 * x;
        } else if (-0.25 <= x && x < 0){
            value = 2.171240593367237661667 + 5.70778203747481756346* x;
        } else if (-0.5 <= x && x < -0.25){
            value = 2.226035845481337351211 + 5.92696304593121632164 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0 <= x && x < 0.25){
            value = - 11.00419101450331474149 ;
        } else if (-0.75 <= x &&  x < -0.5){
            value = - 3.441413161593263778173 ;
        } else if (0.5 <= x && x < 0.75){
            value = - 0.678016042762234026960 ;
        } else if (0.75 <= x && x < 1.0){
            value = 0.096859434680319146709 ;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 ;
        } else if ( 0.25 <= x && x < 0.5){
            value = 2.90038524911627897507 ;
        } else if (-0.25 <= x && x < 0){
            value = 5.70778203747481756346 ;
        } else if (-0.5 <= x && x < -0.25){
            value = 5.92696304593121632164 ;
        } else{
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.5 <= x && x < 1.0) {
            value = 3.464101615137754587055 - 3.464101615137754587055*x;
        } else if(0.0 <= x && x < 0.5){
            value = 3.464101615137754587055 * x ;
        } else {
            value = 0.0;
        }   
    } else if (deriv == 1) {
        if (0.5 <= x && x < 1.0) {
            value =  -3.464101615137754587055;
        } else if(0.0 <= x && x < 0.5){
            value = 3.464101615137754587055;
        } else {
            value = 0.0;
        }   
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_evaluator2(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.25 <= x && x < 0.5) {
            value = 5.89923827193725729450 - 14.13397337635908350175 *x;
        } else if (0.0<= x && x < 0.25){
            value = 9.46297971138994567625 *x;
        } else if (0.75 <= x && x < 1.0) {
            value = -1.677990269774715967092 + 1.677990269774715967092 * x;
        } else if (0.5 <= x && x < 0.75) {
            value = -2.664250113839495385577 + 2.99300339519442185840 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.25 <= x && x < 0.5) {
            value = -14.13397337635908350175;
        } else if (0.0<= x && x < 0.25){
            value = 9.46297971138994567625 ;
        } else if (0.75 <= x && x < 1.0) {
            value =  1.677990269774715967092;
        } else if (0.5 <= x && x < 0.75) {
            value =  2.99300339519442185840 ;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}
    
} // namespace lawa
