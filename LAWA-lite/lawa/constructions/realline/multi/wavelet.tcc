#include <cassert>
#include <iostream>

namespace lawa {
  
template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(int _d)
    : d(_d), vanishingMoments(_d)
{
    _initialize(d);
    if (d > 4) {
        std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
}

template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(const Basis<T,Orthogonal,R,Multi> &_basis)
    : d(_basis.d), vanishingMoments(_basis.d)
{
    _initialize(d);
    if (d > 4) {
        std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet implemented for d = " << d << std::endl;
        exit(1);
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
    //return pow2i<T>(-(j+3));
    return _initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
Wavelet<T,Orthogonal,R,Multi>::getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const
{
    refinement_j = j + _addRefinementLevel;

    long shift = this->_shift(k);
    int  type  = this->_type(k);

    refinement_k_first = _shiftFactor*shift+_offsets[type];
    return &(_refCoeffs[type]);
}

template <typename T>
int
Wavelet<T,Orthogonal,R,Multi>::getRefinementLevel(int j) const
{
    return j + _addRefinementLevel;
}


template <typename T>
void
Wavelet<T,Orthogonal,R,Multi>::_initialize(int d)
{
    switch (d) {
        case 1:
            _numSplines         = 1;
            _initialticsize     = pow2i<T>(-1);
            _addRefinementLevel = 1;
            _shiftFactor        = 2;

            _evaluator = new Evaluator[1];
            _evaluator[0] = _constant_wavelet_inner_evaluator0;

            _support = new Support<T>[1];
            _support[0] = Support<T>(0.,1.);

            _singularSupport = new DenseVector<Array<T> >[1];
            _singularSupport[0].engine().resize(3,0);
            _singularSupport[0] = 0., 0.5, 1.;

            _refCoeffs = new DenseVector<Array<long double> >[1];
            _refCoeffs[0].engine().resize(2,0);
            _refCoeffs[0] =  1.L, -1.L;
            _refCoeffs[0] *= std::pow(2.L,-0.5L);

            _offsets = new long[1];
            _offsets[0] =  0;

            break;

        case 2:
            _numSplines = 3;
            _initialticsize = pow2i<T>(-3);
            _addRefinementLevel = 3;
            _shiftFactor        = 8;

            //Reordering necessary due to old construction is necessary: rhs estimates are for old cons.!
            _evaluator = new Evaluator[3];
            _evaluator[2] = _linear_wavelet_inner_evaluator0;
            _evaluator[1] = _linear_wavelet_inner_evaluator1;
            _evaluator[0] = _linear_wavelet_inner_evaluator2;

            _support = new Support<T>[3];
            _support[2] = Support<T>( 0.,1.);
            _support[1] = Support<T>(-1.,1.);
            _support[0] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[3];
            _singularSupport[2].engine().resize(9,0);
            _singularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[1].engine().resize(15,0);
            _singularSupport[1] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[0].engine().resize(13,0);
            _singularSupport[0] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;

            _max_support = Support<T>(-1,1);

            _refCoeffs = new DenseVector<Array<long double> >[3];
            _refCoeffs[2].engine().resize(7,0);
            _refCoeffs[2] =  0.0704344921276534107756L, -0.576367220716899526482L, 0.639023596671335788701L,
                                  1.81444421027701819162L,   -3.31416570156693953276L,  1.08813353643605850361L,
                                  0.278497086771773164537L;
            _refCoeffs[1].engine().resize(15,0);
            _refCoeffs[1] = -0.0624170162041586273949L, -0.124834032408317254790L,  0.312085081020793136974L,
                                  0.749004194449903528739L,   0.121355651311707938338L, -1.50496515109302569038L,
                                 -0.724655526773638313573L,   0.L,                       2.70859400988290668974L,
                                 -1.42781347611197973525L,   -0.380130236065090770011L,  0.286094144563627797955L,
                                  0.119205893568178249148L,  -0.0476823574272712996592L,-0.0238411787136356498296L;
            _refCoeffs[0].engine().resize(15,0);
            _refCoeffs[0] = -0.0614538064570225674674L, -0.122907612914045134935L, 0.307269032285112837337L,
                                  0.737445677484270809608L,   0.242390522570959039273L,-2.21918643896707489002L,
                                  0.0308223293141810753687L,  2.17124059336723766167L, -1.95533103707150536639L,
                                  0.870385464299548463831L,   0.168829286778896578467L,-0.145289152020478720063L,
                                 -0.0605371466751994666929L,  0.0242148586700797866771L,0.0121074293350398933386L;
            _refCoeffs[0] *= std::pow(2.L,-1.5L);
            _refCoeffs[1] *= std::pow(2.L,-1.5L);
            _refCoeffs[2] *= std::pow(2.L,-1.5L);

            _offsets = new long[3];
            _offsets[2] =  0;
            _offsets[1] = -8;
            _offsets[0] = -8;

            _H1SemiNorms = new long double[3];
            _H1SemiNorms[2] =  std::sqrt(397.584415584415584416L);
            _H1SemiNorms[1] =  std::sqrt(244.857142857142857143L);
            _H1SemiNorms[0] =  std::sqrt(335.558441558441558442L);

            break;

        case 3:
            _numSplines = 6;
            _initialticsize = pow2i<T>(-4);
            _addRefinementLevel = 4;
            _shiftFactor        = 16;

            _evaluator = new Evaluator[6];
            _evaluator[0] = _quadratic_wavelet_inner_evaluator0;
            _evaluator[1] = _quadratic_wavelet_inner_evaluator1;
            _evaluator[2] = _quadratic_wavelet_inner_evaluator2;
            _evaluator[3] = _quadratic_wavelet_inner_evaluator3;
            _evaluator[4] = _quadratic_wavelet_inner_evaluator4;
            _evaluator[5] = _quadratic_wavelet_inner_evaluator5;

            _support = new Support<T>[6];
            _support[0] = Support<T>(0.,1.);
            _support[1] = Support<T>(0.,1.);
            _support[2] = Support<T>(-1.,1.);
            _support[3] = Support<T>(-1.,1.);
            _support[4] = Support<T>(-1.,1.);
            _support[5] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[6];
            _singularSupport[0] = DenseVector<Array<T> >(17);
            _singularSupport[1] = DenseVector<Array<T> >(17);
            _singularSupport[2] = DenseVector<Array<T> >(29);
            _singularSupport[3] = DenseVector<Array<T> >(29);
            _singularSupport[4] = DenseVector<Array<T> >(25);
            _singularSupport[5] = DenseVector<Array<T> >(25);
            _singularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _singularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _singularSupport[2] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _singularSupport[3] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _singularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;

            _max_support = Support<T>(-1,1);

            _refCoeffs = new DenseVector<Array<long double> >[6];
            _refCoeffs[0].engine().resize(14,0);
            _refCoeffs[0] =                            0.0128266866722189260497L, -0.0559640497163821462537L,
                                 0.201215435407263077173L, -2.44304797151868845934L,    5.01689144219717754102L,
                                -1.49983458213893093765L,  -2.67590770526710016897L,    0.486001371590130349617L,
                                 1.35445687638883579514L,  -0.0658794676691616960327L, -0.503186520979130770356L,
                                 0.215040714402484464061L, -0.0453573969046096306311L,  0.00274516753589365617248L;
            _refCoeffs[1].engine().resize(14,0);
            _refCoeffs[1] =                           -0.000915996874565641175157L, 0.0265639127004805607709L,
                                -0.142032263245871410316L,  0.491587510138786841830L,   -1.34292547650285961160L,
                                 1.80303659734814751882L,   0.608131936660204932051L,   -4.61749695265816022830L,
                                 4.54973416715414481731L,  -0.459278023381016642557L,   -1.84270641691420461810L,
                                 1.14698780811210356425L,  -0.227699637379863890094L,    0.00701283484267380709974L;
            _refCoeffs[2].engine().resize(30,0);
            _refCoeffs[2] =                           -0.000587273974134284569257L, -0.00176182192240285370777L,
                                 0.0152691233274913988007L, 0.0505055617755484729561L,   -0.0387727638533256973144L,
                                -0.252565853559131112011L,  0.0109006521064529315779L,    0.751626753143426433452L,
                                 0.131194417340064876211L, -1.75743469498505181457L,     -1.14358520948971039761L,
                                 5.37732404447849841841L,  -2.91552926853981668792L,     -0.619271118755360086863L,
                                 0.L,                       0.L,                          0.724869677753034565042L,
                                -0.273957539933356048382L, -0.424291115783150304622L,     0.612494091159866441960L,
                                -0.134681711544408637645L, -0.105154433065888516400L,    -0.0865358101304027219772L,
                                 0.0277945237247314559225L, 0.0594720551417553510961L,    0.00849678412066896354369L,
                                -0.0125966656855964120689L,-0.00380829427704077574177L,   0.000439418570427781816358L,
                                 0.000146472856809260605453L;
            _refCoeffs[3].engine().resize(30,0);
            _refCoeffs[3] =                              0.0000791262252900211285761L, 0.000237378675870063385728L,
                                -0.00205728185754054934298L, -0.00680485537494181705755L,   0.00523539712290425015158L,
                                 0.0340634756359976522844L,  -0.00176399465305462761871L,  -0.102247013744252589558L,
                                -0.0172365516821472814200L,   0.241734481019712183041L,     0.156253114750098018080L,
                                -0.696643952800272454257L,    0.324713968045916678232L,     0.0644367086364204529505L,
                                 0.L,                         0.L,                          1.01106677157937853260L,
                                -2.01012573985896680233L,    -1.51673479291028773315L,      4.97911911591941212612L,
                                -1.14136144890343037194L,    -1.24434505349016481770L,     -1.01401704218214477329L,
                                 0.325492191000923912041L,    0.696676697419039354043L,     0.0995364770722015527208L,
                                -0.147559417644005964852L,   -0.0446109867295831986763L,    0.00514742154572113830880L,
                                 0.00171580718190704610293L;
            _refCoeffs[4].engine().resize(30,0);
            _refCoeffs[4] =                              0.000102000328969584692298L,  0.000306000986908754076893L,
                                -0.00265200855320920199974L, -0.00877202829138428353759L,   0.00713195296624978726756L,
                                 0.0450599352196930104157L,  -0.0122342185321396614352L,   -0.164750508289248228285L,
                                -0.00816130964527967490745L,  0.501080684146624097214L,     0.152982116115941973510L,
                                 0.484700842065798329957L,   -3.12957282732159598184L,      0.779458241585808449154L,
                                 2.20340313220727341073L,     1.61981639451656894387L,     -3.57997653037485195205L,
                                 0.0948037064793073380219L,   1.34345307579537589384L,     -0.0223214564607085013675L,
                                -0.0265694042646556245245L,  -0.261203868609348449492L,    -0.208022653274724327961L,
                                 0.0666888699496056256548L,   0.142832365015284436264L,     0.0204078319223121038670L,
                                -0.0302515636006807703493L,  -0.00914582155369418638466L,   0.00105528710234932919823L,
                                 0.000351762367449776399410L;
            _refCoeffs[5].engine().resize(30,0);
            _refCoeffs[5] =                              0.0000983648834454673979802L, 0.000295094650336402193940L,
                                -0.00255748696958215234748L, -0.00845937997631019622629L,   0.00671042278164081339608L,
                                 0.0429519213042708765196L,  -0.00744742023225750652765L,  -0.144487601827944335746L,
                                -0.0150946949298144550885L,   0.431884783333983083890L,     0.000655273899011682066927L,
                                 0.445488064897857925478L,   -3.05913557799686498371L,      2.57602406395131441654L,
                                 1.91891257775548164957L,    -4.48546488864417022677L,      3.04637049507008871938L,
                                -0.0368633528973838852982L,  -1.04888016645802703618L,      0.225998341743256388769L,
                                -0.0294486545449911364975L,   0.134278363874709911054L,     0.105736125261204650514L,
                                -0.0338994207331588126176L,  -0.0726025987776971257779,    -0.0103734088724102889667L,
                                 0.0153770862318969752737L,   0.00464888653522466694320L,  -0.000536409984833615416523L,
                                -0.000178803328277871805508L;

            _refCoeffs[0] *= 1./4.;
            _refCoeffs[1] *= 1./4.;
            _refCoeffs[2] *= 1./4.;
            _refCoeffs[3] *= 1./4.;
            _refCoeffs[4] *= 1./4.;
            _refCoeffs[5] *= 1./4.;

            _offsets = new long[6];
            _offsets[0] =  0;
            _offsets[1] =  0;
            _offsets[2] =  -16;
            _offsets[3] =  -16;
            _offsets[4] =  -16;
            _offsets[5] =  -16;

            _H1SemiNorms = new long double[6];
            _H1SemiNorms[0] =  std::sqrt(3770.6381231418184814L/4.L);
            _H1SemiNorms[1] =  std::sqrt(4878.5702463923893863L/4.L);
            _H1SemiNorms[2] =  std::sqrt(964.220175873625172424L);
            _H1SemiNorms[3] =  std::sqrt(3157.160763755615507L/4.L);
            _H1SemiNorms[4] =  std::sqrt(676.111507842511107508L);
            _H1SemiNorms[5] =  std::sqrt(1168.47210336668191741L);

            break;

        case 4:
            _numSplines = 6;
            _initialticsize = pow2i<T>(-4);
            _addRefinementLevel = 3;
            _shiftFactor        = 16;

            _evaluator = new Evaluator[6];
            _evaluator[0] = _cubic_wavelet_inner_evaluator0;
            _evaluator[1] = _cubic_wavelet_inner_evaluator1;
            _evaluator[2] = _cubic_wavelet_inner_evaluator2;
            _evaluator[3] = _cubic_wavelet_inner_evaluator3;
            _evaluator[4] = _cubic_wavelet_inner_evaluator4;
            _evaluator[5] = _cubic_wavelet_inner_evaluator5;

            _support = new Support<T>[6];
            _support[0] = Support<T>(0.,1.);
            _support[1] = Support<T>(0.,1.);
            _support[2] = Support<T>(-1.,1.);
            _support[3] = Support<T>(-1.,1.);
            _support[4] = Support<T>(-1.,1.);
            _support[5] = Support<T>(-1.,1.);

            _singularSupport = new DenseVector<Array<T> >[6];
            _singularSupport[0].engine().resize(9,0);
            _singularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[1].engine().resize(9,0);
            _singularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[2].engine().resize(15,0);
            _singularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[3].engine().resize(15,0);
            _singularSupport[3] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _singularSupport[4].engine().resize(13,0);
            _singularSupport[4] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
            _singularSupport[5].engine().resize(13,0);
            _singularSupport[5] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;

            _max_support = Support<T>(-1,1);

            _refCoeffs = new DenseVector<Array<long double> >[6];
            _refCoeffs[0].engine().resize(14,0);
            _refCoeffs[0] =                            0.211914898968468518860L, -0.0472397067358769825195L,
                                 0.612324176362161983972L, -3.77233294617792048971L,   5.84290173559537136883L,
                                -1.08633452514518673692L,  -3.83459408967896857449L,   1.41141413757794318706L,
                                 1.23008504500417279735L,  -0.0329537785833027628217L,-1.20903898601653731040L,
                                 0.810502565582400589883L, -0.204748785789724792182L,  0.0681002590369992030858L;
            _refCoeffs[1].engine().resize(14,0);
            _refCoeffs[1] =                           -0.0155794646981312991147L, 0.0602797987828925001635L,
                                -0.281572680104588270825L,  0.751940880582288565463L, -1.82766639370529571078L,
                                 1.61758203992707332621L,   2.37943495227619944644L,  -5.76092511515242664165L,
                                 4.26522355075436972153L,   0.746459223273579118041L, -4.60500414433309610986L,
                                 3.34110349656541493881L,  -0.901585828798742936622L,  0.230309684630463352182L;
            _refCoeffs[2].engine().resize(30,0);
            _refCoeffs[2] =                             -0.0195580644521241847758L,-0.0121426891707365422167L,
                                 0.0148307505627752851182L,   0.200956767305319333561L, -0.114474673694457569113L,
                                -0.355474856238127053332L,   -0.281043597782019634878L,  1.46453321678406169808L,
                                -0.000644261226579700550744L,-1.76159554020049779041L,  -2.84870273444244159580L,
                                 6.53096093490833542162L,    -2.73221494717838963219L,  -0.393131334376091888673L,
                                 0.L,                         0.L,                       1.12039566463693547624L,
                                -0.586944670238777824507L,   -2.12549765099280496818L,   2.61121676596332247220L,
                                -0.451642663600581453045L,   -0.146119554450980074420L, -0.389364405226265810917L,
                                 0.159914211093940252717L,    0.157498122840339587243L,  0.0351004964508652265640L,
                                -0.0848810416850084686413L,  -0.00586422489389452393257L,0.00547893342846835013916L,
                                 0.00841104587541561210545L;
            _refCoeffs[3].engine().resize(30,0);
            _refCoeffs[3] =                              0.00826220725232468441722L,0.00511457712212091222684,
                                -0.00629526026040754438076L, -0.0849834387998201084801L, 0.0491865483763372042127L,
                                 0.149897838964449211004L,    0.116439142376403905102L, -0.632763771219819133375L,
                                 0.00146981001896114498526L,  0.773445151779366113379L,  1.22606473850212222355L,
                                -2.70987276244970375721L,     0.994855942890322344397L,  0.109179275447342800168L,
                                 0.L,                         0.L,                       1.68355902093375845067L,
                                -1.69047766435055296490L,    -4.52937869119672571389L,   6.35667059024853567043L,
                                -0.955974255435916400604L,   -0.520875091193768100325L, -1.16318849147368099258L,
                                 0.475112706234268715587L,    0.468576565283153625306L,  0.104716089970779217169L,
                                -0.252608244390480100688L,   -0.0174598203231073028214L, 0.0162986877660911226172L,
                                 0.0250285979276447740279L;
            _refCoeffs[4].engine().resize(30,0);
            _refCoeffs[4] =                              0.00118662293125793290344L,0.000532475676933906769009L,
                                -0.00130829450864805226886L, -0.0134178669760076200341L, 0.0181787550189524438421L,
                                 0.0178909227036158521662L,  -0.0139935316066808033859L,-0.280027409295067248865L,
                                 0.00735978672926965722201L,  0.540088046374024923158L,  0.282730538954526642976L,
                                 0.901164634807985088283L,   -3.59597672997964773665L,   0.367714234348273582742L,
                                 2.42024378272409704942L,     1.79272690850629438757L,  -3.11569858296215764266L,
                                -1.02494353303302156504L,     2.21431696671276166169L,  -0.373601505934363244108L,
                                 0.468515034866185402063L,   -0.445669100075347689932L, -0.542575652269063608370L,
                                 0.216318163121035402632L,    0.214650890310242641257L,  0.0485550820576930980042L,
                                -0.115873453384063683873L,   -0.00802468851782075487594L,0.00746257921991479996764L,
                                 0.0114749234788251774056L;
            _refCoeffs[5].engine().resize(30,0);
            _refCoeffs[5] =                              0.00166325905412781376306L, 0.000924939932864417174325L,
                                -0.00147663824252679317747L, -0.0177359871606028896369L,  0.0156586658575696059346L,
                                 0.0282917321704542054793L,   0.00753014546516630945249L,-0.225354162346751353400L,
                                -0.0155707144779878903312L,   0.453673559221570769471L,  -0.0530005881765670224092L,
                                 1.07509212556543765569L,    -4.05351836832630334422L,    2.24340518699856692127L,
                                 4.28385475732783247240L,    -6.10504394267743267937L,    2.59500728001207241285L,
                                 1.15399825985332605319L,    -2.16658707631630021693L,    0.838199856026823196218L,
                                -0.427013769645531489758L,    0.281500473707012175040L,   0.280116173948380893546L,
                                -0.111841690620988501923L,   -0.110938480161408450900L,  -0.0250765089995902168078L,
                                 0.0598822517026479662611L,   0.00414658810391109897501L,-0.00385701687190878417950L,
                                -0.00593031092386433366700L;

            _refCoeffs[0] *= std::pow(2.L,-1.5L);
            _refCoeffs[1] *= std::pow(2.L,-1.5L);
            _refCoeffs[2] *= std::pow(2.L,-1.5L);
            _refCoeffs[3] *= std::pow(2.L,-1.5L);
            _refCoeffs[4] *= std::pow(2.L,-1.5L);
            _refCoeffs[5] *= std::pow(2.L,-1.5L);

            _offsets = new long[6];
            _offsets[0] =  0;
            _offsets[1] =  0;
            _offsets[2] =  -16;
            _offsets[3] =  -16;
            _offsets[4] =  -16;
            _offsets[5] =  -16;

            _H1SemiNorms = new long double[6];
            _H1SemiNorms[0] =  std::sqrt(3581.1787249033985863L/4.L);
            _H1SemiNorms[1] =  std::sqrt(4745.4284412264144244L/4.L);
            _H1SemiNorms[2] =  std::sqrt(3551.9683102854048364L/4.L);
            _H1SemiNorms[3] =  std::sqrt(3470.3742489363095669L/4.L);
            _H1SemiNorms[4] =  std::sqrt(2382.5255297084867072L/4.L);
            _H1SemiNorms[5] =  std::sqrt(4381.1340895992035858L/4.L);

            break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
    return;
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
    return k>=0 ? k % _numSplines : _numSplines - (_numSplines-1-k)% _numSplines - 1;
}

} // namespace lawa
