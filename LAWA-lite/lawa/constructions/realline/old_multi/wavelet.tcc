#include <cassert>
#include <iostream>

namespace lawa {
  
template <typename T>
    T
    _linear_wavelet_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_evaluator2(T x, unsigned short deriv);

//------------------------------------------------------------------------------

template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(int _d)
    : d(_d), vanishingMoments(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 2: _numSplines = 3;

                _evaluator = new Evaluator[3];
                _evaluator[0] = _linear_wavelet_evaluator0;
                _evaluator[1] = _linear_wavelet_evaluator1;                
                _evaluator[2] = _linear_wavelet_evaluator2;
            
                _support = new Support<T>[3];
                _support[0] = Support<T>(-1,1);
                _support[1] = Support<T>(-1,1);
                _support[2] = Support<T>( 0,1);
            
                _singularSupport = new DenseVector<Array<T> >[3];
                _singularSupport[0].engine().resize(13,0);
                _singularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[1].engine().resize(13,0);
                _singularSupport[1] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[2] = linspace(0.0,1.0,9);

                _max_support = Support<T>(-1,1);
                break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}

template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(const Basis<T,Orthogonal,R,Multi> &_basis)
    : d(_basis.d), vanishingMoments(_basis.d)
{
    assert(d>=2);

    switch (d) {
        case 2: _numSplines = 3;

                _evaluator = new Evaluator[3];
                _evaluator[0] = _linear_wavelet_evaluator0;
                _evaluator[1] = _linear_wavelet_evaluator1;
                _evaluator[2] = _linear_wavelet_evaluator2;

                _support = new Support<T>[3];
                _support[0] = Support<T>(-1,1);
                _support[1] = Support<T>(-1,1);
                _support[2] = Support<T>( 0,1);

                _singularSupport = new DenseVector<Array<T> >[3];
                _singularSupport[0].engine().resize(13,0);
                _singularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[1].engine().resize(13,0);
                _singularSupport[1] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[2] = linspace(0.0,1.0,9);

                _max_support = Support<T>(-1,1);
                break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,R,Multi>::~Wavelet()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
Wavelet<T,Orthogonal,R,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2ih<T>(2*j*deriv+j) *
    _evaluator[type](pow2i<T>(j)*x - shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,R,Multi>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
Support<T>
Wavelet<T,Orthogonal,R,Multi>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,R,Multi>::singularSupport(int j, long k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;
}

template <typename T>
T
Wavelet<T,Orthogonal,R,Multi>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

template <typename T>
long
Wavelet<T,Orthogonal,R,Multi>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
Wavelet<T,Orthogonal,R,Multi>::_type(long k) const
{
    return k>=0 ? (int) (k%3) : (int) _numSplines - (int)((-k+2)% ((int)_numSplines)) - 1;
}

//------------------------------------------------------------------------------
    
template <typename T>
T
_linear_wavelet_evaluator0(T _x, unsigned short deriv)
{
    long double value = 0.0L;
    long double x = (long double) _x;
    if (deriv == 0) {
        if (0L <= x && x < 0.125L) {
            value = 2.17124059336723766167L - 33.0125730435099442245L * x;
        } else if(-0.375L <= x && x < -0.25L) {
            value = -7.14234036204314274860L - 19.69261569230427143433L * x;
        } else if(0.25L <= x && x < 0.375L){
            value = 2.27349781934085223456L - 5.6124494201652150829L * x;
        } else if(-0.5L <= x && x < -0.375L){
            value = -1.24277494216897627173L - 3.96044123930649416268L * x;
        } else if(0.375L <= x && x < 0.5L ){
            value = 1.11118460317702247406L - 2.51294751039500238824L * x;
        } else if(-1.0L <= x && x < -0.75L){
            value = -0.491630451656180539739L - 0.491630451656180539739L * x;
        } else if(0.75L <= x && x < 1L){
            value = 0.096859434680319146709L - 0.096859434680319146709L * x;
        } else if(0.5L <= x && x < 0.75L) {
            value = -0.48429717340159573354L + 0.67801604276223402696L * x;
        } else if(-0.75L <= x && x < -0.5L){
            value = 2.458152258280902698695L + 3.44141316159326377817L *x;
        } else if(-0.125L <= x && x < 0L){
            value = 2.17124059336723766167L + 17.1233461124244526904L *x;
        } else if(-0.25L <= x && x < -0.125L){
            value = 2.28083109759543704076L + 18.0000701462500477231L *x;
        } else if(0.125L <= x && x < 0.25L){
            value = -4.78104753844255919661L + 22.6057320109684306418L* x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0L <= x && x < 0.125L) {
            value = - 33.0125730435099442245L ;
        } else if(-0.375L <= x && x < -0.25L) {
            value = - 19.69261569230427143433L ;
        } else if(0.25L <= x && x < 0.375L){
            value = -5.6124494201652150829L ;
        } else if(-0.5L <= x && x < -0.375L){
            value = -3.96044123930649416268L ;
        } else if(0.375L <= x && x < 0.5L ){
            value = -2.51294751039500238824L ;
        } else if(-1.0L <= x && x < -0.75L){
            value = -0.491630451656180539739L ;
        } else if(0.75L <= x && x < 1L){
            value = -0.096859434680319146709L ;
        } else if(0.5L <= x && x < 0.75L) {
            value =  0.67801604276223402696L ;
        } else if(-0.75L <= x && x <-0.5L){
            value =  3.44141316159326377817L ;
        } else if(-0.125L <= x && x < 0L){
            value = 17.1233461124244526904L ;
        } else if(-0.25L <= x && x < -0.125L){
            value = 18.0000701462500477231L ;
        } else if(0.125L <= x && x < 0.25L){
            value = 22.6057320109684306418L;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return (T)value;
}

template <typename T>
T
_linear_wavelet_evaluator1(T _x, unsigned short deriv)
{
    long double value = 0.0L;
    long double x = (long double) _x;
    if (deriv == 0) {
        if (0.125L <= x && x< 0.25L) {
            value = 6.84500149587779311472L - 33.09125988795909139990L* x;
        } else if(-0.375L <= x && x < -0.25L){
            value = -4.757606755902492947821L - 13.01056641923786902976L *x;
        } else if(-0.5L <= x && x < -0.375L){
            value = -1.761589978102878832866L - 5.02118834510556472321L*x;
        } else if(0.5L <= x && x < 0.75L){
            value = 0.953647148545425993185L - 1.335106007963596390459L*x;
        } else if(-1.0L <= x && x < -0.75L){
            value = -0.4993361296332690191592L - 0.4993361296332690191592L*x;
        } else if(0.75L <= x && x < 1.0L){
            value = -0.1907294297090851986370L + 0.1907294297090851986370L*x;
        } else if(-0.125L <= x && x < 0.0L){
            value = 5.797244214189106508582L*x;
        } else if(0.0 <= x && x < 0.125){
            value = 21.66875207906325351789L* x;
        } else if(-0.75L <= x && x < -0.5L){
            value = 2.496680648166345095796L + 3.495352907432883134114L*x;
        } else if(0.375L <= x && x < 0.5L){
            value = -2.378803377951246473911L + 5.32979504502974854373L*x;
        } else if(-0.25L <= x && x < -0.125L){
            value = 0.0556540975457490632366L + 6.24247699455509901447L*x;
        } else if(0.25L <= x && x < 0.375L){
            value = -3.523179956205757665732L + 8.38146592037511172192L*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.125L <= x && x< 0.25L) {
            value = - 33.09125988795909139990L;
        } else if(-0.375L <= x && x < -0.25L){
            value =  - 13.01056641923786902976L ;
        } else if(-0.5L <= x && x < -0.375L){
            value =  - 5.02118834510556472321L;
        } else if(0.5L <= x && x < 0.75L){
            value = - 1.335106007963596390459L;
        } else if(-1.0L <= x && x < -0.75L){
            value = - 0.4993361296332690191592L;
        } else if(0.75L <= x && x < 1.0L){
            value = 0.1907294297090851986370L;
        } else if(-0.125L <= x && x <= 0.0L){//TODO: //ttttttttttest <=0.0
            value = 5.797244214189106508582L;
        } else if(0.0L <= x && x < 0.125L){
            value = 21.66875207906325351789L;
        } else if(-0.75L <= x && x < -0.5L){
            value =  3.495352907432883134114L;
        } else if(0.375L <= x && x < 0.5L){
            value =  5.32979504502974854373L;
        } else if(-0.25L <= x && x < -0.125L){
            value =  6.24247699455509901447L;
        } else if(0.25L <= x && x < 0.375L){
            value =  8.38146592037511172192L;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return (T)value;
}

template <typename T>
T
_linear_wavelet_evaluator2(T _x, unsigned short deriv)
{
    long double value = 0.0L;
    long double x = (long double) _x;
    if (deriv == 0) {
        if (0.625L <= x && x < 0.75L) {
            value = 25.3256618915819297146L - 35.2183939040239842910L*x;
        } else if(0.25L <= x && x < 0.375L){
            value = 3.007148855493370156849L - 9.72312653910588252146L*x;
        } else if(0.375L <= x && x < 0.5L){
            value = 2.887238244145711420071L - 9.40336490884545922339L*x;
        } else if(0.0L <= x && x < 0.125L){
            value = -0.56347593702122728620L*x;
        } else if (0.875L <= x && x < 1.0L){
            value = -2.22797669417418531629L + 2.22797669417418531629L *x;
        } else if (0.125L <= x && x < 0.25L){
            value = -0.717236204972206348034L + 5.17441370275642349806L*x;
        } else if (0.75L <= x && x < 0.875L){
            value = -5.94595223442177053803L + 6.47709159731428271256L*x;
        } else if (0.5L <= x && x < 0.625L){
            value = -22.32888385765284908917L + 41.0288792947516617951L*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.625L <= x && x < 0.75L) {
            value =  - 35.2183939040239842910L;
        } else if(0.25L <= x && x < 0.375L){
            value = - 9.72312653910588252146L;
        } else if(0.375L <= x && x < 0.5L){
            value = - 9.40336490884545922339L;
        } else if(0.0L <= x && x < 0.125L){
            value = -0.56347593702122728620L;
        } else if (0.875L <= x && x < 1.0L){
            value = 2.22797669417418531629L ;
        } else if (0.125L <= x && x < 0.25L){
            value = 5.17441370275642349806L;
        } else if (0.75L <= x && x < 0.875L){
            value = 6.47709159731428271256L;
        } else if (0.5L <= x && x < 0.625L){
            value = 41.0288792947516617951L;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return (T)value;
}

} // namespace lawa
