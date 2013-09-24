#include <cassert>
#include <iostream>

namespace lawa {
    
template <typename T>
BSpline<T,Orthogonal,R,Multi>::BSpline(const int _d)
    : d(_d)
{
    if (d > 4) {
        std::cerr << "BSpline<T,Orthogonal,R,Multi> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
    
    switch (d) {
        case 1:
            _initialticsize     = 1.;
            _numSplines         = 1;
            _addRefinementLevel = 1;
            _shiftFactor        = 2;

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
            _initialticsize     = pow2i<T>(-2);
            _numSplines         = 3;
            _addRefinementLevel = 2;
            _shiftFactor        = 4;
            

            // Order of evaluators interchanged to remain consistend with older lawa versions!!
            _evaluator = new Evaluator[3];
            _evaluator[1] = _linear_bspline_inner_evaluator0;
            _evaluator[2] = _linear_bspline_inner_evaluator1;
            _evaluator[0] = _linear_bspline_inner_evaluator2;

            _support = new Support<T>[3];
            _support[1] = Support<T>( 0.,1.);
            _support[2] = Support<T>( 0.,1);
            _support[0] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[3];
            _singularSupport[1].engine().resize(3,0);
            _singularSupport[1] = 0., 0.5, 1.;
            _singularSupport[2].engine().resize(5,0);
            _singularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[0].engine().resize(9,0);
            _singularSupport[0] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;

            _max_support = Support<T>(-1,1);



            _refCoeffs = new DenseVector<Array<long double> >[3];
            _refCoeffs[1].engine().resize(3,0);
            _refCoeffs[1] = std::sqrt(3.L)/2.L, std::sqrt(3.L), std::sqrt(3.L)/2.L;
            _refCoeffs[2].engine().resize(3,0);
            _refCoeffs[2] = 2.36574492784748641906L, -1.16774841624228445637L, -0.419497567443678991773L;
            _refCoeffs[0].engine().resize(7,0);
            _refCoeffs[0] = 0.122907612914045134935L,-0.737445677484270809608L, 0.744295083998533270801L,
                             2.17124059336723766167L, -0.579807160258591023705L, 0.145289152020478720063L,
                            -0.0242148586700797866771L;
            _refCoeffs[1] *= 0.5L;
            _refCoeffs[2] *= 0.5L;
            _refCoeffs[0] *= 0.5;

            _offsets = new long[3];
            _offsets[1] =  0;
            _offsets[2] =  0;
            _offsets[0] =  -4;

            _H1SemiNorms = new long double[3];
            _H1SemiNorms[1] = std::sqrt(12.L);
            _H1SemiNorms[2] = std::sqrt(75.2727272727272727273L);;
            _H1SemiNorms[0] = std::sqrt(52.4415584415584415584L);

            break;

        case 3:
            _initialticsize = pow2i<T>(-3);
            _numSplines = 6;
            _addRefinementLevel = 3;
            _shiftFactor        = 8;

            _evaluator = new Evaluator[6];
            _evaluator[0] = _quadratic_bspline_inner_evaluator0;
            _evaluator[1] = _quadratic_bspline_inner_evaluator1;
            _evaluator[2] = _quadratic_bspline_inner_evaluator2;
            _evaluator[3] = _quadratic_bspline_inner_evaluator3;
            _evaluator[4] = _quadratic_bspline_inner_evaluator4;
            _evaluator[5] = _quadratic_bspline_inner_evaluator5;

            _support = new Support<T>[6];
            _support[0] = Support<T>( 0.,1.);
            _support[1] = Support<T>( 0.,1.);
            _support[2] = Support<T>( 0.,1.);
            _support[3] = Support<T>( 0.,1.);
            _support[4] = Support<T>(-1.,1.);
            _support[5] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[6];
            _singularSupport[0].engine().resize(5,0);
            _singularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[1].engine().resize(5,0);
            _singularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[2].engine().resize(9,0);
            _singularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[3].engine().resize(9,0);
            _singularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[4].engine().resize(17,0);
            _singularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[5].engine().resize(17,0);
            _singularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;

            _max_support = Support<T>(-1,1);

            _refCoeffs = new DenseVector<Array<long double> >[6];
            _refCoeffs[0].engine().resize(6,0);
            _refCoeffs[0] =                                    std::sqrt(15.L/23.L)/2.L, 3.L*std::sqrt(15.L/23.L)/2.L,
                                 2.L*std::sqrt(15.L/23.L), 2.L*std::sqrt(15.L/23.L),     3.L*std::sqrt(15.L/23.L)/2.L,
                                     std::sqrt(15.L/23.L)/2.L;
            _refCoeffs[1].engine().resize(6,0);
            _refCoeffs[1] =                                    std::sqrt(3.L/2.L)/2.L, 3.L*std::sqrt(3.L/2.L)/2.L,
                                  std::sqrt(3.L/2.L),         -std::sqrt(3.L/2.L),    -3.L*std::sqrt(3.L/2.L)/2.L,
                                  -std::sqrt(3.L/2.L)/2.L;
            _refCoeffs[2].engine().resize(6,0);
            _refCoeffs[2] =                                -1.32957575266269475104L, 1.00928819006999249208L,
                                 1.11372017505655343851L,  -3.57158506454512406189L, 2.64084214163580360767L,
                                 0.614683388767643909165L;
            _refCoeffs[3].engine().resize(6,0);
            _refCoeffs[3] =                                 3.35387951165860404196L, -0.141229992246845435862L,
                                -1.29954452099952139115L,  -0.599442481500795701774L, 1.08068229875787216245L,
                                 0.312647416473029867432L;
            _refCoeffs[4].engine().resize(14,0);
            _refCoeffs[4] =                                -0.000661624087242333507038L, 0.0191870985300276717041L,
                                -0.103410901756548045803L,  0.378884765460252640419L,   -1.07547200857749368479L,
                                 1.10969884455047509484L,   2.08619476135110666967L,     2.08619476135110666967L,
                                -0.795461007681955947706L,  0.000218300557375081050850L, 0.247205601265849432686L,
                                -0.146044023363594721918L,  0.0292055710131498807034L,  -0.00100708865562585795529L;
            _refCoeffs[5].engine().resize(14,0);
            _refCoeffs[5] =                                -0.00113465969390614975730L,  0.0329051311232783429617L,
                                -0.176288182431553100420L,  0.619104902646686873282L,   -1.70854321322203183450L,
                                 2.15119201578740835235L,   1.45657142326475405947L,    -3.40474081654023230539L,
                                 1.23029213497100169954L,  -0.000708279550556974276803L,-0.379632133141489874567L,
                                 0.224287719297057272008L, -0.0448526868857883735848L,   0.00154664437537201288224L;

            _refCoeffs[0] *= std::pow(2.L,-1.5L);
            _refCoeffs[1] *= std::pow(2.L,-1.5L);
            _refCoeffs[2] *= std::pow(2.L,-1.5L);
            _refCoeffs[3] *= std::pow(2.L,-1.5L);
            _refCoeffs[4] *= std::pow(2.L,-1.5L);
            _refCoeffs[5] *= std::pow(2.L,-1.5L);

            _offsets = new long[6];
            _offsets[0] =  0;
            _offsets[1] =  0;
            _offsets[2] =  0;
            _offsets[3] =  0;
            _offsets[4] =  -8;
            _offsets[5] =  -8;

            _H1SemiNorms = new long double[6];
            _H1SemiNorms[0] =  std::sqrt(13.9130434782608695652L);
            _H1SemiNorms[1] =  8.L;
            _H1SemiNorms[2] =  std::sqrt(268.675020609247651090L);
            _H1SemiNorms[3] =  std::sqrt(131.345247650197830689L);
            _H1SemiNorms[4] =  std::sqrt(81.2723006501163910509L);
            _H1SemiNorms[5] =  std::sqrt(263.708111955787835019L);

            break;
            
        case 4:
            // inner part
            _initialticsize     = pow2i<T>(-2);
            _numSplines         = 6;
            _addRefinementLevel = 2;
            _shiftFactor        = 8;

            _evaluator = new Evaluator[6];
            _evaluator[0] = _cubic_bspline_inner_evaluator0;
            _evaluator[1] = _cubic_bspline_inner_evaluator1;
            _evaluator[2] = _cubic_bspline_inner_evaluator2;
            _evaluator[3] = _cubic_bspline_inner_evaluator3;
            _evaluator[4] = _cubic_bspline_inner_evaluator4;
            _evaluator[5] = _cubic_bspline_inner_evaluator5;

            _support = new Support<T>[6];
            _support[0] = Support<T>(0.,1.);
            _support[1] = Support<T>(0.,1.);
            _support[2] = Support<T>(0.,1.);
            _support[3] = Support<T>(0.,1.);
            _support[4] = Support<T>(-1.,1.);
            _support[5] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[6];
            _singularSupport[0].engine().resize(3,0);
            _singularSupport[0] = 0., 0.5, 1.;
            _singularSupport[1].engine().resize(3,0);
            _singularSupport[1] = 0., 0.5, 1.;
            _singularSupport[2].engine().resize(5,0);
            _singularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[3].engine().resize(5,0);
            _singularSupport[3] = 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[4].engine().resize(9,0);
            _singularSupport[4] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            _singularSupport[5].engine().resize(9,0);
            _singularSupport[5] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;

            _max_support = Support<T>(-1,1);

            _refCoeffs = new DenseVector<Array<long double> >[6];
            _refCoeffs[0].engine().resize(6,0);
            _refCoeffs[0] =                                    std::sqrt(35.L/13.L)/4.L, 3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L),       std::sqrt(35.L/13.L),     3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L)/4.L;
            _refCoeffs[1].engine().resize(6,0);
            _refCoeffs[1] =                                     std::sqrt(35.L/3.L)/4.L,   std::sqrt(35.L/3.L)/2.L,
                                   std::sqrt(35.L/3.L)/2.L,    -std::sqrt(35.L/3.L)/2.L,  -std::sqrt(35.L/3.L)/2.L,
                                  -std::sqrt(35.L/3.L)/4.L;
            _refCoeffs[2].engine().resize(6,0);
            _refCoeffs[2] =                                -1.20968608645812591492L,   0.723769207408211462783L,
                                 2.78393856347267693448L,  -5.12943299583678560134L,   3.02163500540779053256L,
                                 0.583314706422121720912L;
            _refCoeffs[3].engine().resize(6,0);
            _refCoeffs[3] =                                 3.23422528793847648562L,   0.295496503819005929826L,
                                 -2.64548042111713192471L,  0.336105140833428426303L,  0.872145478151140075300L,
                                  0.546845047769619692603L;
            _refCoeffs[4].engine().resize(14,0);
            _refCoeffs[4] =                                 -0.00757430944378488466295L,  0.0304795243208619817595L,
                                 -0.172103404134940755815L,  0.560461227654867482281L,   -1.42578946715790175367L,
                                  0.970460974994275945818L,  2.18881605627051402000L,     2.18881605627051402000L,
                                 -0.619734802216407440776L, -0.418091821163284038735L,    0.735614181183121958188L,
                                 -0.475237518862456267198L,  0.135438721938371823978L,   -0.0366114708090683515524L;
            _refCoeffs[5].engine().resize(14,0);
            _refCoeffs[5] =                                 -0.0183274146720283930429L,   0.0714158787346878889722L,
                                 -0.346389724181854150965L,  0.968551554773324101920L,   -2.38163361520215938629L,
                                  1.98214605626068896928L,   3.24419184855680683262L,    -4.62339059465912210039L,
                                  1.03257150206081408287L,   0.612554526218050167521L,   -1.11025696825143994205L,
                                  0.717898961851500034340L, -0.204656443589658444638L,    0.0553244321003903398448L;

            _refCoeffs[0] *= 0.5L;
            _refCoeffs[1] *= 0.5L;
            _refCoeffs[2] *= 0.5L;
            _refCoeffs[3] *= 0.5L;
            _refCoeffs[4] *= 0.5L;
            _refCoeffs[5] *= 0.5L;

            _offsets = new long[6];
            _offsets[0] =  0;
            _offsets[1] =  0;
            _offsets[2] =  0;
            _offsets[3] =  0;
            _offsets[4] = -8;
            _offsets[5] = -8;

            _H1SemiNorms = new long double[6];
            _H1SemiNorms[0] =  std::sqrt(12.9230769230769230769L);
            _H1SemiNorms[1] =  std::sqrt(56.L);
            _H1SemiNorms[2] =  std::sqrt(259.228164412374246759L);
            _H1SemiNorms[3] =  std::sqrt(123.272155615509401277L);
            _H1SemiNorms[4] =  std::sqrt(79.7761595473175714506L);
            _H1SemiNorms[5] =  std::sqrt(258.536491525265347505L);

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
    delete[] _refCoeffs;
    delete[] _offsets;
    delete[] _H1SemiNorms;
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
    //return pow2i<T>(-(j+3));
    return _initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
BSpline<T,Orthogonal,R,Multi>::getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const
{
    refinement_j = j + _addRefinementLevel;

    long shift = this->_shift(k);
    int  type  = this->_type(k);

    refinement_k_first = _shiftFactor*shift+_offsets[type];
    return &(_refCoeffs[type]);
}

template <typename T>
int
BSpline<T,Orthogonal,R,Multi>::getRefinementLevel(int j) const
{
    return j + _addRefinementLevel;
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
    //return k>=0 ? (int) k%3 : (int) _numSplines - (int)(-k+2)%((int)_numSplines) - 1;
    return k>=0 ? k % _numSplines : _numSplines - (_numSplines-1-k) % _numSplines - 1;

}
    
} // namespace lawa
