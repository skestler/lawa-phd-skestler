#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC 1

#include <cassert>
#include <lawa/math/linspace.h>
#include <lawa/math/pow2.h>
#include <lawa/constructions/interval/multi/_constant_evaluator.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Orthogonal,Interval,Multi>::Basis(int _d, int j)
    : mra(_d, j), d(_d), d_(_d), j0(mra.j0), _bc(2,0), _j(j0), psi(*this), refinementbasis(d)
{
    assert(d>=1);
    
    /* L_2 orthonormal multiwavelet bases without Dirichlet boundary conditions are not
     * implemented yet. They require __different__ boundary adapted wavelets and scaling functions.
     */

    _numLeftParts = 0;
    _numRightParts = 0;

    if (d==1) {
        _addRefinementLevel = 1;
        _shiftFactor        = 2;
        //left part
        _numLeftParts = 1;
        _leftEvaluator = new Evaluator[1];
        _leftEvaluator[0] = _constant_wavelet_inner_evaluator0;

        _leftSupport = new Support<T>[1];
        _leftSupport[0] = Support<T>(0.,1.);

        _leftSingularSupport = new DenseVector<Array<T> >[1];
        _leftSingularSupport[0].engine().resize(3,0);
        _leftSingularSupport[0] = 0., 0.5, 1.;

        _leftRefCoeffs = new DenseVector<Array<long double> >[1];
        _leftRefCoeffs[0].engine().resize(2,0);
        _leftRefCoeffs[0] =  1.L, -1.L;

        _leftRefCoeffs[0] *= std::pow(2.L,-0.5L); // valued multiplied by (-1)... difference to mathematica?

        _leftOffsets = new long[1];
        _leftOffsets[0] =  0;

        //_leftH1SemiNorms = new long double[0];

        //inner part
        _numInnerParts = 1;
        _innerEvaluator = new Evaluator[1];
        _innerEvaluator[0] = _constant_wavelet_inner_evaluator0;

        _innerSupport = new Support<T>[1];
        _innerSupport[0] = Support<T>(0.,1.);

        _innerSingularSupport = new DenseVector<Array<T> >[1];
        _innerSingularSupport[0].engine().resize(3,0);
        _innerSingularSupport[0] = 0., 0.5, 1.;

        _innerRefCoeffs = new DenseVector<Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize(2,0);
        _innerRefCoeffs[0] =  1.L, -1.L;
        _innerRefCoeffs[0] *= std::pow(2.L,-0.5L);

        _innerOffsets = new long[1];
        _innerOffsets[0] =  0;

        //_innerH1SemiNorms = new long double[0];

        //right part
        _numRightParts = 0;
        //_rightEvaluator = new Evaluator[0];
        //_rightSupport = new Support<T>[0];
        //_rightSingularSupport = new DenseVector<Array<T> >[0];
        //_rightRefCoeffs = new DenseVector<Array<long double> >[0];
        //_rightOffsets = new long[0];
        //_rightH1SemiNorms = new long double[0];
    }
    else {
        this->enforceBoundaryCondition<DirichletBC>();
    }

    setLevel(_j);
}

template <typename T>
Basis<T,Orthogonal,Interval,Multi>::~Basis()
{
    if (d==1) {
        delete[] _leftEvaluator;
        delete[] _innerEvaluator;
        delete[] _leftSupport;
        delete[] _innerSupport;
        delete[] _leftSingularSupport;
        delete[] _innerSingularSupport;
        delete[] _leftRefCoeffs;
        delete[] _innerRefCoeffs;
        delete[] _leftOffsets;
        delete[] _innerOffsets;
    }
    else {
        if (_numLeftParts>0) {
            delete[] _leftEvaluator;
            delete[] _leftSupport;
            delete[] _leftSingularSupport;
            delete[] _leftRefCoeffs;
            delete[] _leftOffsets;
            delete[] _leftH1SemiNorms;
        }
        if (_numInnerParts>0) {
            delete[] _innerEvaluator;
            delete[] _innerSupport;
            delete[] _innerSingularSupport;
            delete[] _innerRefCoeffs;
            delete[] _innerOffsets;
            delete[] _innerH1SemiNorms;
        }
        if (_numRightParts>0) {
            delete[] _rightEvaluator;
            delete[] _rightSupport;
            delete[] _rightSingularSupport;
            delete[] _rightRefCoeffs;
            delete[] _rightOffsets;
            delete[] _rightH1SemiNorms;
        }
    }
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    _bc = 1,1;

    switch (d) {
        case 2:
            _addRefinementLevel = 3;
            _shiftFactor        = 8;
            //left part
            _numLeftParts = 2;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_wavelet_left_evaluator0;
            _leftEvaluator[1] = _linear_wavelet_inner_evaluator0;

            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);

            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(9,0);
            _leftSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;

            _leftRefCoeffs = new DenseVector<Array<long double> >[2];
            _leftRefCoeffs[0].engine().resize(7,0);
            _leftRefCoeffs[0] =  3.33475836931829630962L,  -1.75789096554035096452L, -0.468007564634539555940L,
                                 0.352232501259116508276L,  0.146763542191298545115L,-0.0587054168765194180459L,
                                -0.0293527084382597090230L;
            _leftRefCoeffs[1].engine().resize(7,0);
            _leftRefCoeffs[1] = 0.070434492127653410776L, -0.57636722071689952648L, 0.63902359667133578870L,
                                1.8144442102770181916L,   -3.3141657015669395328L,  1.0881335364360585036L,
                                0.27849708677177316454L;

            _leftRefCoeffs[0] *= -std::pow(2.L,-1.5L); // valued multiplied by (-1)... difference to mathematica?
            _leftRefCoeffs[1] *=  std::pow(2.L,-1.5L);

            _leftOffsets = new long[2];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;

            _leftH1SemiNorms = new long double[2];
            _leftH1SemiNorms[0] =  std::sqrt(315.827533859762919478L);
            _leftH1SemiNorms[1] =  std::sqrt(397.584415584415584416L);;

            //inner part
            _numInnerParts = 3;
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _linear_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _linear_wavelet_inner_evaluator2;

            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(-1.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);

            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0].engine().resize(9,0);
            _innerSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[1].engine().resize(15,0);
            _innerSingularSupport[1] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[2].engine().resize(13,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;


            _innerRefCoeffs = new DenseVector<Array<long double> >[3];
            _innerRefCoeffs[0].engine().resize(7,0);
            _innerRefCoeffs[0] =  0.0704344921276534107756L, -0.576367220716899526482L, 0.639023596671335788701L,
                                  1.81444421027701819162L,   -3.31416570156693953276L,  1.08813353643605850361L,
                                  0.278497086771773164537L;
            _innerRefCoeffs[1].engine().resize(15,0);
            _innerRefCoeffs[1] = -0.0624170162041586273949L, -0.124834032408317254790L,  0.312085081020793136974L,
                                  0.749004194449903528739L,   0.121355651311707938338L, -1.50496515109302569038L,
                                 -0.724655526773638313573L,   0.L,                       2.70859400988290668974L,
                                 -1.42781347611197973525L,   -0.380130236065090770011L,  0.286094144563627797955L,
                                  0.119205893568178249148L,  -0.0476823574272712996592L,-0.0238411787136356498296L;
            _innerRefCoeffs[2].engine().resize(15,0);
            _innerRefCoeffs[2] = -0.0614538064570225674674L, -0.122907612914045134935L, 0.307269032285112837337L,
                                  0.737445677484270809608L,   0.242390522570959039273L,-2.21918643896707489002L,
                                  0.0308223293141810753687L,  2.17124059336723766167L, -1.95533103707150536639L,
                                  0.870385464299548463831L,   0.168829286778896578467L,-0.145289152020478720063L,
                                 -0.0605371466751994666929L,  0.0242148586700797866771L,0.0121074293350398933386L;
            _innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[2] *= std::pow(2.L,-1.5L);

            _innerOffsets = new long[3];
            _innerOffsets[0] =  0;
            _innerOffsets[1] = -8;
            _innerOffsets[2] = -8;

            _innerH1SemiNorms = new long double[3];
            _innerH1SemiNorms[0] =  std::sqrt(397.584415584415584416L);
            _innerH1SemiNorms[1] =  std::sqrt(244.857142857142857143L);
            _innerH1SemiNorms[2] =  std::sqrt(335.558441558441558442L);

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _linear_wavelet_right_evaluator0;

            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);

            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(7,0);
            _rightRefCoeffs[0] = -0.107000114803417747921L, -0.214000229606835495841L, 0.535000574017088739604L,
                                  1.28400137764101297505L,   0.208037317578970289244L,-2.57992857933775636329L,
                                 -1.24226099344595652887L;
            _rightRefCoeffs[0] *= -std::pow(2.L,-1.5L);

            _rightOffsets = new long[1];
            _rightOffsets[0] =  0;

            _rightH1SemiNorms = new long double[1];
            _rightH1SemiNorms[0] =  std::sqrt(107.263375231146171431L);
            break;
        
        case 3:
            _addRefinementLevel = 4;
            _shiftFactor        = 16;
            //left part
            _numLeftParts = 4;
            _leftEvaluator = new Evaluator[4];
            _leftEvaluator[0] = _quadratic_wavelet_left_evaluator0;
            _leftEvaluator[1] = _quadratic_wavelet_left_evaluator1;
            _leftEvaluator[2] = _quadratic_wavelet_inner_evaluator0;
            _leftEvaluator[3] = _quadratic_wavelet_inner_evaluator1;

            _leftSupport = new Support<T>[4];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);

            _leftSingularSupport = new DenseVector<Array<T> >[4];
            _leftSingularSupport[0] = DenseVector<Array<T> >(13);
            _leftSingularSupport[1] = DenseVector<Array<T> >(17);
            _leftSingularSupport[2] = DenseVector<Array<T> >(17);
            _leftSingularSupport[3] = DenseVector<Array<T> >(17);
            _leftSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[2] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[3] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;

            _leftRefCoeffs = new DenseVector<Array<long double> >[4];
            _leftRefCoeffs[0].engine().resize(15,0);
            _leftRefCoeffs[0] = 7.34142886366690571014L,  -4.32324127529253091787L,   -0.0148993032084588272842L,
                                1.42933391151209638388L,  -0.234378132475796052981L,   0.0230059106395800543011L,
                               -0.204317957585265375598L, -0.161347143157051860204L,   0.0517405162918780304874L,
                                0.110799856937372102120L,  0.0158308787794303546941L, -0.0234673317737972813416L,
                               -0.00709477472231080598698L,0.000818627852574323767729L,0.000272875950858107922576L;
            _leftRefCoeffs[1].engine().resize(15,0);
            _leftRefCoeffs[1] = -0.175226007747828497337L,  1.33231359644221527055L,   -2.03213893910277911952L,
                                -1.64437451962213360457L,   5.01683485527306782513L,   -1.14725876761424721952L,
                                -1.23061673230483248656L,  -1.00325061275701909267L,    0.322043497348364334898L,
                                 0.689287250020311826693L,  0.0984806452588233827096L, -0.145994383728320622226L,
                                -0.0441378369411201881149L, 0.00509282733936002170556L, 0.00169760911312000723519;
            _leftRefCoeffs[2].engine().resize(14,0);
            _leftRefCoeffs[2] =                             0.0128266866722189260497L, -0.0559640497163821462537L,
                                 0.201215435407263077173L, -2.44304797151868845934L,    5.01689144219717754102L,
                                -1.49983458213893093765L,  -2.67590770526710016897L,    0.486001371590130349617L,
                                 1.35445687638883579514L,  -0.0658794676691616960327L, -0.503186520979130770356L,
                                 0.215040714402484464061L, -0.0453573969046096306311L,  0.00274516753589365617248L;
            _leftRefCoeffs[3].engine().resize(14,0);
            _leftRefCoeffs[3] =                            -0.000915996874565641175157L, 0.0265639127004805607709L,
                                -0.142032263245871410316L,  0.491587510138786841830L,   -1.34292547650285961160L,
                                 1.80303659734814751882L,   0.608131936660204932051L,   -4.61749695265816022830L,
                                 4.54973416715414481731L,  -0.459278023381016642557L,   -1.84270641691420461810L,
                                 1.14698780811210356425L,  -0.227699637379863890094L,    0.00701283484267380709974L;

            _leftRefCoeffs[0] *= 1./4.;
            _leftRefCoeffs[1] *= 1./4.;
            _leftRefCoeffs[2] *= 1./4.;
            _leftRefCoeffs[3] *= 1./4.;

            _leftOffsets = new long[4];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;

            _leftH1SemiNorms = new long double[4];
            _leftH1SemiNorms[0] =  std::sqrt(1689.1831469636316641L);
            _leftH1SemiNorms[1] =  std::sqrt(815.437862271242019281L);
            _leftH1SemiNorms[2] =  std::sqrt(3770.6381231418184814L/4.L);
            _leftH1SemiNorms[3] =  std::sqrt(4878.5702463923893863L/4.L);


            //inner part
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _quadratic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _quadratic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _quadratic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _quadratic_wavelet_inner_evaluator3;
            _innerEvaluator[4] = _quadratic_wavelet_inner_evaluator4;
            _innerEvaluator[5] = _quadratic_wavelet_inner_evaluator5;

            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            _innerSupport[3] = Support<T>(-1.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);

            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0] = DenseVector<Array<T> >(17);
            _innerSingularSupport[1] = DenseVector<Array<T> >(17);
            _innerSingularSupport[2] = DenseVector<Array<T> >(29);
            _innerSingularSupport[3] = DenseVector<Array<T> >(29);
            _innerSingularSupport[4] = DenseVector<Array<T> >(25);
            _innerSingularSupport[5] = DenseVector<Array<T> >(25);
            _innerSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[2] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[3] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;

            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(14,0);
            _innerRefCoeffs[0] =                            0.0128266866722189260497L, -0.0559640497163821462537L,
                                 0.201215435407263077173L, -2.44304797151868845934L,    5.01689144219717754102L,
                                -1.49983458213893093765L,  -2.67590770526710016897L,    0.486001371590130349617L,
                                 1.35445687638883579514L,  -0.0658794676691616960327L, -0.503186520979130770356L,
                                 0.215040714402484464061L, -0.0453573969046096306311L,  0.00274516753589365617248L;
            _innerRefCoeffs[1].engine().resize(14,0);
            _innerRefCoeffs[1] =                           -0.000915996874565641175157L, 0.0265639127004805607709L,
                                -0.142032263245871410316L,  0.491587510138786841830L,   -1.34292547650285961160L,
                                 1.80303659734814751882L,   0.608131936660204932051L,   -4.61749695265816022830L,
                                 4.54973416715414481731L,  -0.459278023381016642557L,   -1.84270641691420461810L,
                                 1.14698780811210356425L,  -0.227699637379863890094L,    0.00701283484267380709974L;
            _innerRefCoeffs[2].engine().resize(30,0);
            _innerRefCoeffs[2] =                           -0.000587273974134284569257L, -0.00176182192240285370777L,
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
            _innerRefCoeffs[3].engine().resize(30,0);
            _innerRefCoeffs[3] =                              0.0000791262252900211285761L, 0.000237378675870063385728L,
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
            _innerRefCoeffs[4].engine().resize(30,0);
            _innerRefCoeffs[4] =                              0.000102000328969584692298L,  0.000306000986908754076893L,
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
            _innerRefCoeffs[5].engine().resize(30,0);
            _innerRefCoeffs[5] =                              0.0000983648834454673979802L, 0.000295094650336402193940L,
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

            _innerRefCoeffs[0] *= 1./4.;
            _innerRefCoeffs[1] *= 1./4.;
            _innerRefCoeffs[2] *= 1./4.;
            _innerRefCoeffs[3] *= 1./4.;
            _innerRefCoeffs[4] *= 1./4.;
            _innerRefCoeffs[5] *= 1./4.;

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  -15;
            _innerOffsets[3] =  -15;
            _innerOffsets[4] =  -15;
            _innerOffsets[5] =  -15;

            _innerH1SemiNorms = new long double[6];
            _innerH1SemiNorms[0] =  std::sqrt(3770.6381231418184814L/4.L);
            _innerH1SemiNorms[1] =  std::sqrt(4878.5702463923893863L/4.L);
            _innerH1SemiNorms[2] =  std::sqrt(964.220175873625172424L);
            _innerH1SemiNorms[3] =  std::sqrt(3157.160763755615507L/4.L);
            _innerH1SemiNorms[4] =  std::sqrt(676.111507842511107508L);
            _innerH1SemiNorms[5] =  std::sqrt(1168.47210336668191741L);


            //right part
            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _quadratic_wavelet_right_evaluator0;
            _rightEvaluator[1] = _quadratic_wavelet_right_evaluator1;

            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(0.,1.);
            _rightSupport[1] = Support<T>(0.,1.);

            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0] = DenseVector<Array<T> >(13);
            _rightSingularSupport[1] = DenseVector<Array<T> >(13);
            _rightSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _rightSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[2];
            _rightRefCoeffs[0].engine().resize(15,0);
            _rightRefCoeffs[0] =                              -0.000114497443892427926522L, -0.000343492331677283779566L,
                                  0.00297693354120312608957L,  0.00984678017474880168089L,  -0.00790337721930663062058,
                                 -0.0502735386409631708149L,   0.0110711479726297769007L,    0.176130682621472212526L,
                                  0.0141356360586135612533L,  -0.547135481139192757034L,     0.00430442950944302670782L,
                                 -0.834704980616085529054L,    4.37995567073963987782L,     -3.02869064594082812035L,
                                 -3.17209055617314295964L;
            _rightRefCoeffs[1].engine().resize(15,0);
            _rightRefCoeffs[1] =                              -0.000582268855038003371693L, -0.00174680656511401011508L,
                                  0.0151389902309880876640L,   0.0500751215332682899656L,   -0.0384058632272029463584L,
                                 -0.250303964050425621308L,    0.00985992805274688018376L,   0.742085813082314558117L,
                                  0.131553702458883543274L,   -1.72847682305559402948L,     -1.11668290900582672902L,
                                  5.38485064046866963707L,    -3.00904975376343648823L,     -0.796018915521131065822L,
                                 -0.110216370589382144526L;

            _rightRefCoeffs[0] *= 1./4.;
            _rightRefCoeffs[1] *= 1./4.;

            _rightOffsets = new long[2];
            _rightOffsets[0] =  1;
            _rightOffsets[1] =  1;

            _rightH1SemiNorms = new long double[2];
            _rightH1SemiNorms[0] =  std::sqrt(872.85470590556468096L);
            _rightH1SemiNorms[1] =  std::sqrt(950.23153980918785971L);

            break;
            
        case 4:
            _addRefinementLevel = 3;
            _shiftFactor        = 16;
            // left part
            _numLeftParts = 4;
            _leftEvaluator = new Evaluator[4];
            _leftEvaluator[0] = _cubic_wavelet_left_evaluator0;
            _leftEvaluator[1] = _cubic_wavelet_left_evaluator1;
            _leftEvaluator[2] = _cubic_wavelet_inner_evaluator0;
            _leftEvaluator[3] = _cubic_wavelet_inner_evaluator1;
            
            _leftSupport = new Support<T>[4];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[4];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(9,0);
            _leftSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[2].engine().resize(9,0);
            _leftSingularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[3].engine().resize(9,0);
            _leftSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[4];
            _leftRefCoeffs[0].engine().resize(15,0);
            _leftRefCoeffs[0] =  9.96825463930175006856L,  -3.64662794301701212754L,  -1.72146472727434846315L,
                                 2.96847209211156307389L,  -1.07608624861786460327,    0.578622693490629381008L,
                                -0.397264462426789392875L, -0.400953974066674688364L,  0.160462581368940836059L,
                                 0.159072235220550385378L,  0.0359146177147943384943L,-0.0858526536425712577082,
                                -0.00594379904285631537199L,0.00553075754282213126303L,0.00850265706425028894903L;
            _leftRefCoeffs[1].engine().resize(15,0);
            _leftRefCoeffs[1] = -0.227957578161734593738L,  2.30332425068216343077L,  -1.68253494167435384206L,
                                -5.20257457174657321107L,   6.88037200400476413447L,  -1.09982321724826908893,
                                -0.488736025233599518193L, -1.16197270348703346833L,   0.475517036761813332104L,
                                 0.468752917383957175504L,  0.104655980664234953224L, -0.252676836677631112457L,
                                -0.0174618882848239483496L, 0.0163054523354829342512,  0.0250363964778949084260L;
            _leftRefCoeffs[2].engine().resize(14,0);
            _leftRefCoeffs[2] =                             0.211914898968468518860L, -0.0472397067358769825195L,
                                 0.612324176362161983972L, -3.77233294617792048971L,   5.84290173559537136883L,
                                -1.08633452514518673692L,  -3.83459408967896857449L,   1.41141413757794318706L,
                                 1.23008504500417279735L,  -0.0329537785833027628217L,-1.20903898601653731040L,
                                 0.810502565582400589883L, -0.204748785789724792182L,  0.0681002590369992030858L;
            _leftRefCoeffs[3].engine().resize(14,0);
            _leftRefCoeffs[3] =                            -0.0155794646981312991147L, 0.0602797987828925001635L,
                                -0.281572680104588270825L,  0.751940880582288565463L, -1.82766639370529571078L,
                                 1.61758203992707332621L,   2.37943495227619944644L,  -5.76092511515242664165L,
                                 4.26522355075436972153L,   0.746459223273579118041L, -4.60500414433309610986L,
                                 3.34110349656541493881L,  -0.901585828798742936622L,  0.230309684630463352182L;

            _leftRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[2] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[3] *= std::pow(2.L,-1.5L);

            _leftOffsets = new long[4];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;

            _leftH1SemiNorms = new long double[4];
            _leftH1SemiNorms[0] =  std::sqrt(5472.168230418815658L/4.L);
            _leftH1SemiNorms[1] =  std::sqrt(3658.4430353007925789L/4.L);
            _leftH1SemiNorms[2] =  std::sqrt(3581.1787249033985872L/4.L);
            _leftH1SemiNorms[3] =  std::sqrt(4745.4284412264144235L/4.L);

            // inner part
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _cubic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _cubic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _cubic_wavelet_inner_evaluator3;
            _innerEvaluator[4] = _cubic_wavelet_inner_evaluator4;
            _innerEvaluator[5] = _cubic_wavelet_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            _innerSupport[3] = Support<T>(-1.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(9,0);
            _innerSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[1].engine().resize(9,0);
            _innerSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[2].engine().resize(15,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[3].engine().resize(15,0);
            _innerSingularSupport[3] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[4].engine().resize(13,0);
            _innerSingularSupport[4] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
            _innerSingularSupport[5].engine().resize(13,0);
            _innerSingularSupport[5] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(14,0);
            _innerRefCoeffs[0] =                            0.211914898968468518860L, -0.0472397067358769825195L,
                                 0.612324176362161983972L, -3.77233294617792048971L,   5.84290173559537136883L,
                                -1.08633452514518673692L,  -3.83459408967896857449L,   1.41141413757794318706L,
                                 1.23008504500417279735L,  -0.0329537785833027628217L,-1.20903898601653731040L,
                                 0.810502565582400589883L, -0.204748785789724792182L,  0.0681002590369992030858L;
            _innerRefCoeffs[1].engine().resize(14,0);
            _innerRefCoeffs[1] =                           -0.0155794646981312991147L, 0.0602797987828925001635L,
                                -0.281572680104588270825L,  0.751940880582288565463L, -1.82766639370529571078L,
                                 1.61758203992707332621L,   2.37943495227619944644L,  -5.76092511515242664165L,
                                 4.26522355075436972153L,   0.746459223273579118041L, -4.60500414433309610986L,
                                 3.34110349656541493881L,  -0.901585828798742936622L,  0.230309684630463352182L;
            _innerRefCoeffs[2].engine().resize(30,0);
            _innerRefCoeffs[2] =                             -0.0195580644521241847758L,-0.0121426891707365422167L,
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
            _innerRefCoeffs[3].engine().resize(30,0);
            _innerRefCoeffs[3] =                              0.00826220725232468441722L,0.00511457712212091222684,
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
            _innerRefCoeffs[4].engine().resize(30,0);
            _innerRefCoeffs[4] =                              0.00118662293125793290344L,0.000532475676933906769009L,
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
            _innerRefCoeffs[5].engine().resize(30,0);
            _innerRefCoeffs[5] =                              0.00166325905412781376306L, 0.000924939932864417174325L,
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

            _innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[2] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[3] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[4] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[5] *= std::pow(2.L,-1.5L);

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  -15;
            _innerOffsets[3] =  -15;
            _innerOffsets[4] =  -15;
            _innerOffsets[5] =  -15;

            _innerH1SemiNorms = new long double[6];
            _innerH1SemiNorms[0] =  std::sqrt(3581.1787249033985863L/4.L);
            _innerH1SemiNorms[1] =  std::sqrt(4745.4284412264144244L/4.L);
            _innerH1SemiNorms[2] =  std::sqrt(3551.9683102854048364L/4.L);
            _innerH1SemiNorms[3] =  std::sqrt(3470.3742489363095669L/4.L);
            _innerH1SemiNorms[4] =  std::sqrt(2382.5255297084867072L/4.L);
            _innerH1SemiNorms[5] =  std::sqrt(4381.1340895992035858L/4.L);

            //right part
            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _cubic_wavelet_right_evaluator0;
            _rightEvaluator[1] = _cubic_wavelet_right_evaluator1;
            
            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(0.,1.);
            _rightSupport[1] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;
            _rightSingularSupport[1].engine().resize(7,0);
            _rightSingularSupport[1] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[2];
            _rightRefCoeffs[0].engine().resize(15,0);
            _rightRefCoeffs[0] =                              -0.00175688408219078023460L, -0.000925754508436546215605L,
                                  0.00166225914750846803799L,  0.0190418501000516459908L,  -0.0193588653236211675258L,
                                 -0.0289617697271932935231L,  -0.000163958707092606003783L, 0.286009738020544576093L,
                                  0.0169323716197137782115L,  -0.597167034247085329268L,    0.0783897634859167994357L,
                                 -1.61759750023053098540L,     5.62030470810278152907L,    -2.73895672808327636225L,
                                 -5.72183982291327789481L;
            _rightRefCoeffs[1].engine().resize(15,0);
            _rightRefCoeffs[1] =                              -0.0210717289670313478080,   -0.0130965475120609487646,
                                  0.0159503629099407980868,    0.216424920466073176724,    -0.122559045517901730267,
                                 -0.383239957986316512419,    -0.304936904470756387579,     1.56468573108351611838,
                                  0.000994248161337121536308, -1.87272034507059529094,     -3.03700555827392189530,
                                  7.05112499271747390693,     -3.01562383975439626345,     -0.539402266004325008587,
                                 -0.107513637558473482228;

            _rightRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _rightRefCoeffs[1] *= std::pow(2.L,-1.5L);

            _rightOffsets = new long[2];
            _rightOffsets[0] =  1;
            _rightOffsets[1] =  1;

            _rightH1SemiNorms = new long double[2];
            _rightH1SemiNorms[0] =  std::sqrt(3577.8921753706540727L/4.L);
            _rightH1SemiNorms[1] =  std::sqrt(3494.6338698754857046L/4.L);
            break;
            
        default: std::cerr << "BC for Wavelet<T,Orthogonal,Interval,Multi> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    refinementbasis.enforceBoundaryCondition<BC>();
    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Orthogonal,Interval,Multi> &
Basis<T,Orthogonal,Interval,Multi>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJ(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<int>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJI(int j) const
{
    assert(j>=j0);
    return (pow2i<int>(j)-1)*_numInnerParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJ(int j) const
{
    assert(j>=j0);
    return Range<int>(1,cardJ(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<int>(1,cardJL());
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJI(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+1,cardJL()+cardJI(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+cardJI(j)+1,cardJ(j));
}

template <typename T>
void
Basis<T,Orthogonal,Interval,Multi>::getLowerEnclosingScaling(int j_wavelet, long k_wavelet,
                                                             int &j_scaling, long &k_scaling) const
{
    j_scaling = j_wavelet;
    if (k_wavelet<=rangeJL(j_wavelet).lastIndex()) {
        k_scaling = mra.rangeI(j_scaling).firstIndex();
        return;
    }
    if (k_wavelet>=rangeJR(j_wavelet).firstIndex()) {
        k_scaling = mra.rangeI(j_scaling).lastIndex();
        return;
    }
    long k_tilde = (k_wavelet / _numInnerParts + 1) * mra._numInnerParts;
    long k_scaling_first = k_tilde - 2* mra._numInnerParts - 2;
    k_scaling_first = std::max(k_scaling_first, (long)mra.rangeI(j_scaling).firstIndex());
    long k_scaling_last  = k_tilde + 2* mra._numInnerParts + 2;
    k_scaling_last  = std::min(k_scaling_last, (long)mra.rangeI(j_scaling).lastIndex());
    k_scaling = (k_scaling_last + k_scaling_first) / 2;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,Multi>::getLowerEnclosingWavelet(int j_wavelet1, long k_wavelet1,
                                                             int &j_wavelet2, long &k_wavelet2) const
{
    j_wavelet2 = j_wavelet1 - 1;
    if (k_wavelet1<=this->rangeJL(j_wavelet1).lastIndex()) {
        k_wavelet2 = 1;
        return;
    }
    if (k_wavelet1>=this->rangeJR(j_wavelet1).firstIndex()) {
        k_wavelet2 = rangeJR(j_wavelet2).lastIndex();
        return;
    }
    long k_tilde = k_wavelet1 / 2;
    long k_wavelet_first = k_tilde - 2*_numInnerParts;
    k_wavelet_first = std::max(k_wavelet_first, (long)rangeJ(j_wavelet2).firstIndex());
    long k_wavelet_last  = k_tilde + 2*_numInnerParts;
    k_wavelet_last  = std::min(k_wavelet_last,  (long)rangeJ(j_wavelet2).lastIndex());
    k_wavelet2 = (k_wavelet_last + k_wavelet_first) / 2;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getScalingNeighborsForScaling(int j_scaling1, long k_scaling1,
                                const SecondBasis &secondbasis,
                                int &j_scaling2, long &k_scaling_first, long &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_scaling2 = j_scaling1;
    k_scaling_first = k_scaling1 - 2* mra._numInnerParts;
    k_scaling_first = std::max(k_scaling_first, (long)mra.rangeI(j_scaling2).firstIndex());
    k_scaling_last  = k_scaling1 + 2* mra._numInnerParts;
    k_scaling_last  = std::min(k_scaling_last, (long)mra.rangeI(j_scaling2).lastIndex());
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getWaveletNeighborsForScaling(int j_scaling, long k_scaling,
                                const SecondBasis &secondbasis,
                                int &j_wavelet, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_wavelet = j_scaling;
    k_wavelet_first = k_scaling - 2*_numInnerParts;
    k_wavelet_first = std::max(k_wavelet_first, (long)rangeJ(j_wavelet).firstIndex());
    k_wavelet_last  = k_scaling + 2*_numInnerParts;
    k_wavelet_last  = std::min(k_wavelet_last, (long)rangeJ(j_wavelet).lastIndex());
    return;
}

template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getBSplineNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                const SecondRefinementBasis &secondrefinementbasis,
                                int &j_bspline, long &k_bspline_first, long &k_bspline_last) const
{
    ct_assert(SecondRefinementBasis::Side==Orthogonal and SecondRefinementBasis::Domain==Interval
              and SecondRefinementBasis::Cons==MultiRefinement);
    j_bspline = j_wavelet + _addRefinementLevel - 1;
    Support<T> supp = psi.support(j_wavelet,k_wavelet);
    T h = 1.L/(T)(refinementbasis.mra.cardI(j_bspline)-1);
    k_bspline_first = std::floor(supp.l1/h) - refinementbasis.mra._numLeftParts;
    k_bspline_first = std::max(k_bspline_first,(long)refinementbasis.mra.rangeI(j_bspline).firstIndex());
    k_bspline_last  = std::ceil(supp.l2/h) + refinementbasis.mra._numRightParts;
    k_bspline_last = std::min(k_bspline_last,(long)refinementbasis.mra.rangeI(j_bspline).lastIndex());
    //std::cerr << "(" << j_wavelet << "," << k_wavelet << "): " << std::endl;
    //std::cerr << "   "  << supp << " " << refinementbasis.mra.phi.support(j_bspline,k_bspline_first)
    //          << " " << refinementbasis.mra.phi.support(j_bspline,k_bspline_last) << std::endl;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getScalingNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                const SecondBasis &secondbasis,
                                int &j_scaling, long &k_scaling_first, long &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_scaling = j_wavelet;
    long k_tilde = (k_wavelet / _numInnerParts + 1) * mra._numInnerParts;
    k_scaling_first = k_tilde - 2* mra._numInnerParts - 2;
    k_scaling_first = std::max(k_scaling_first, (long)mra.rangeI(j_scaling).firstIndex());
    k_scaling_last  = k_tilde + 2* mra._numInnerParts + 2;
    k_scaling_last  = std::min(k_scaling_last, (long)mra.rangeI(j_scaling).lastIndex());
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1, const SecondBasis &secondbasis,
                                int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1;
    long k_tilde = k_wavelet1;
    k_wavelet_first = k_tilde - 2*_numInnerParts;
    k_wavelet_first = std::max(k_wavelet_first, (long)rangeJ(j_wavelet2).firstIndex());
    k_wavelet_last  = k_tilde + 2*_numInnerParts;
    k_wavelet_last  = std::min(k_wavelet_last,  (long)rangeJ(j_wavelet2).lastIndex());
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getLowerWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                     const SecondBasis &secondbasis,
                                     int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1-1;
    long k_tilde = k_wavelet1 / 2;
    k_wavelet_first = k_tilde - 2*_numInnerParts;
    k_wavelet_first = std::max(k_wavelet_first, (long)rangeJ(j_wavelet2).firstIndex());
    k_wavelet_last  = k_tilde + 2*_numInnerParts;
    k_wavelet_last  = std::min(k_wavelet_last,  (long)rangeJ(j_wavelet2).lastIndex());
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,Multi>
::getHigherWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                     const SecondBasis &secondbasis,
                                     int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1+1;
    long k_tilde = 2*k_wavelet1;
    k_wavelet_first = k_tilde - 4*_numInnerParts;
    k_wavelet_first = std::max(k_wavelet_first, (long)rangeJ(j_wavelet2).firstIndex());
    k_wavelet_last  = k_tilde + 3*_numInnerParts;
    k_wavelet_last  = std::min(k_wavelet_last,  (long)rangeJ(j_wavelet2).lastIndex());
    return;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC
