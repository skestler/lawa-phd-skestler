#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_MRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_MRA_TCC 1

#include <cassert>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::MRA(int _d, int j)
    : d(_d), j0((j==-1) ? 0 : j), _bc(2,0), _j(j0), phi(*this)
{
    assert(d>=1);
    assert(j0>=0);
    
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
        _leftEvaluator[0] = _constant_bspline_inner_evaluator0;

        _leftSupport = new Support<T>[1];
        _leftSupport[0] = Support<T>(0.,1.);

        _leftSingularSupport = new DenseVector<Array<T> >[1];
        _leftSingularSupport[0].engine().resize(2,0);
        _leftSingularSupport[0] = 0., 1.;

        _leftRefCoeffs = new DenseVector<Array<long double> >[1];
        _leftRefCoeffs[0].engine().resize(2,0);
        _leftRefCoeffs[0] = 1.L, 1.L;

        _leftRefCoeffs[0] *= std::pow(2.L,-0.5L);
        _leftRefCoeffs[1] *= std::pow(2.L,-0.5L);

        _leftOffsets = new long[1];
        _leftOffsets[0] =  0;

        //_leftH1SemiNorms = new long double[0];

        //inner part
        _numInnerParts = 1;
        _innerEvaluator = new Evaluator[1];
        _innerEvaluator[0] = _constant_bspline_inner_evaluator0;

        _innerSupport = new Support<T>[1];
        _innerSupport[0] = Support<T>(0.,1.);

        _innerSingularSupport = new DenseVector<Array<T> >[1];
        _innerSingularSupport[0].engine().resize(2,0);
        _innerSingularSupport[0] = 0., 1.;

        _innerRefCoeffs = new DenseVector<Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize(2,0);
        _innerRefCoeffs[0] = 1.L, 1.L;
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
        //_rightOffsets   = new long[0];
        //_rightH1SemiNorms = new long double[0];
    }
    else {
        this->enforceBoundaryCondition<DirichletBC>();
    }

    setLevel(_j);
}

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::~MRA()
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

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardI(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<int>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardIL(int /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardII(int j) const
{
    assert(j>=j0);
    return (pow2i<int>(j)-1)*_numInnerParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardIR(int /*j*/) const
{
    return _numRightParts;
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeI(int j) const
{
    assert(j>=j0);
    return Range<int>(0,cardI(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIL(int /*j*/) const
{
    return Range<int>(0,cardIL() - 1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeII(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL(), cardIL()+cardII(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL()+cardII(j),cardI(j)-1);
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

    switch (d) {
        case 2: 
            _addRefinementLevel = 2;
            _shiftFactor        = 4;
            //left part
            _numLeftParts = 2;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_bspline_inner_evaluator0;
            _leftEvaluator[1] = _linear_bspline_inner_evaluator1;
            
            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0].engine().resize(3,0);
            _leftSingularSupport[0] = 0., 0.5, 1.;
            _leftSingularSupport[1].engine().resize(5,0);
            _leftSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[2];
            _leftRefCoeffs[0].engine().resize(3,0);
            _leftRefCoeffs[0] = std::sqrt(3.L)/2.L, std::sqrt(3.L), std::sqrt(3.L)/2.L;
            _leftRefCoeffs[1].engine().resize(3,0);
            _leftRefCoeffs[1] = 2.36574492784748641906L, -1.16774841624228445637L, -0.419497567443678991773L;

            _leftRefCoeffs[0] *= 0.5L;
            _leftRefCoeffs[1] *= 0.5L;

            _leftOffsets = new long[2];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;

            _leftH1SemiNorms = new long double[2];
            _leftH1SemiNorms[0] = std::sqrt(12.L);
            _leftH1SemiNorms[1] = std::sqrt(75.2727272727272727273L);

            //inner part
            _numInnerParts = 3;
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_bspline_inner_evaluator0;
            _innerEvaluator[1] = _linear_bspline_inner_evaluator1;
            _innerEvaluator[2] = _linear_bspline_inner_evaluator2;
            
            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 0.5, 1.;
            _innerSingularSupport[1].engine().resize(5,0);
            _innerSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[2].engine().resize(9,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[3];
            _innerRefCoeffs[0].engine().resize(3,0);
            _innerRefCoeffs[0] = std::sqrt(3.L)/2.L, std::sqrt(3.L), std::sqrt(3.L)/2.L;
            _innerRefCoeffs[1].engine().resize(3,0);
            _innerRefCoeffs[1] = 2.36574492784748641906L, -1.16774841624228445637L, -0.419497567443678991773L;
            _innerRefCoeffs[2].engine().resize(7,0);
            _innerRefCoeffs[2] = 0.122907612914045134935L,-0.737445677484270809608L, 0.744295083998533270801L,
                                 2.17124059336723766167L, -0.579807160258591023705L, 0.145289152020478720063L,
                                -0.0242148586700797866771L;
            _innerRefCoeffs[0] *= 0.5L;
            _innerRefCoeffs[1] *= 0.5L;
            _innerRefCoeffs[2] *= 0.5;

            _innerOffsets = new long[3];
            _innerOffsets[0] =  0;
            _innerOffsets[1] =  0;
            _innerOffsets[2] = -4;

            _innerH1SemiNorms = new long double[3];
            _innerH1SemiNorms[0] = std::sqrt(12.L);
            _innerH1SemiNorms[1] = std::sqrt(75.2727272727272727273L);
            _innerH1SemiNorms[2] = std::sqrt(52.4415584415584415584L);


            //right part
            _numRightParts = 0;
            //_rightEvaluator = new Evaluator[0];
            //_rightSupport = new Support<T>[0];
            //_rightSingularSupport = new DenseVector<Array<T> >[0];

            //_rightRefCoeffs = new DenseVector<Array<long double> >[0];
            //_rightOffsets   = new long[0];
            //_rightH1SemiNorms = new long double[0];
            break;
            
        case 3:
            _addRefinementLevel = 3;
            _shiftFactor        = 8;
            // left part
            _numLeftParts = 5;
            _leftEvaluator = new Evaluator[5];
            _leftEvaluator[0] = _quadratic_bspline_left_evaluator0;
            _leftEvaluator[1] = _quadratic_bspline_inner_evaluator0;
            _leftEvaluator[2] = _quadratic_bspline_inner_evaluator1;
            _leftEvaluator[3] = _quadratic_bspline_inner_evaluator2;
            _leftEvaluator[4] = _quadratic_bspline_inner_evaluator3;
            
            _leftSupport = new Support<T>[5];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            _leftSupport[4] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[5];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(5,0);
            _leftSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[2].engine().resize(5,0);
            _leftSingularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[3].engine().resize(9,0);
            _leftSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[4].engine().resize(9,0);
            _leftSingularSupport[4] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[5];
            _leftRefCoeffs[0].engine().resize(7,0);
            _leftRefCoeffs[0] = 5.55117912386048439353L,  -1.96151750303782862164L,  0.00138479559372387609709L,
                                0.603401146175499957922L, -0.356497187379049387413L, 0.0712917932693664232410L,
                               -0.00245833769894366976693L;

            _leftRefCoeffs[1].engine().resize(6,0);
            _leftRefCoeffs[1] =                               std::sqrt(15.L/23.L)/2.L, 3.L*std::sqrt(15.L/23.L)/2.L,
                                2.L*std::sqrt(15.L/23.L), 2.L*std::sqrt(15.L/23.L),     3.L*std::sqrt(15.L/23.L)/2.L,
                                std::sqrt(15.L/23.L)/2.L;
            _leftRefCoeffs[2].engine().resize(6,0);
            _leftRefCoeffs[2] =                               std::sqrt(3.L/2.L)/2.L, 3.L*std::sqrt(3.L/2.L)/2.L,
                                std::sqrt(3.L/2.L),         -std::sqrt(3.L/2.L),    -3.L*std::sqrt(3.L/2.L)/2.L,
                               -std::sqrt(3.L/2.L)/2.L;
            _leftRefCoeffs[3].engine().resize(6,0);
            _leftRefCoeffs[3] =                           -1.32957575266269475104L, 1.00928819006999249208L,
                                1.11372017505655343851L,  -3.57158506454512406189L, 2.64084214163580360767L,
                                0.614683388767643909165L;
            _leftRefCoeffs[4].engine().resize(6,0);
            _leftRefCoeffs[4] =                            3.35387951165860404196L, -0.141229992246845435862L,
                               -1.29954452099952139115L,  -0.599442481500795701774L, 1.08068229875787216245L,
                                0.312647416473029867432L;

            _leftRefCoeffs[0] *= std::pow(2.L,-1.5L);;
            _leftRefCoeffs[1] *= std::pow(2.L,-1.5L);;
            _leftRefCoeffs[2] *= std::pow(2.L,-1.5L);;
            _leftRefCoeffs[3] *= std::pow(2.L,-1.5L);;
            _leftRefCoeffs[4] *= std::pow(2.L,-1.5L);;

            _leftOffsets = new long[5];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  1;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;
            _leftOffsets[4] =  1;

            _leftH1SemiNorms = new long double[5];
            _leftH1SemiNorms[0] =  std::sqrt(396.800602767473266279L);
            _leftH1SemiNorms[1] =  std::sqrt(13.9130434782608695652L);
            _leftH1SemiNorms[2] =  8.L;
            _leftH1SemiNorms[3] =  std::sqrt(268.675020609247651090L);
            _leftH1SemiNorms[4] =  std::sqrt(131.345247650197830689L);

            // inner part
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _quadratic_bspline_inner_evaluator0;
            _innerEvaluator[1] = _quadratic_bspline_inner_evaluator1;
            _innerEvaluator[2] = _quadratic_bspline_inner_evaluator2;
            _innerEvaluator[3] = _quadratic_bspline_inner_evaluator3;
            _innerEvaluator[4] = _quadratic_bspline_inner_evaluator4;
            _innerEvaluator[5] = _quadratic_bspline_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(0.,1.);
            _innerSupport[3] = Support<T>(0.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(5,0);
            _innerSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[1].engine().resize(5,0);
            _innerSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[2].engine().resize(9,0);
            _innerSingularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[3].engine().resize(9,0);
            _innerSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[4].engine().resize(17,0);
            _innerSingularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[5].engine().resize(17,0);
            _innerSingularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(6,0);
            _innerRefCoeffs[0] =                               std::sqrt(15.L/23.L)/2.L, 3.L*std::sqrt(15.L/23.L)/2.L,
                                 2.L*std::sqrt(15.L/23.L), 2.L*std::sqrt(15.L/23.L),     3.L*std::sqrt(15.L/23.L)/2.L,
                                     std::sqrt(15.L/23.L)/2.L;
            _innerRefCoeffs[1].engine().resize(6,0);
            _innerRefCoeffs[1] =                               std::sqrt(3.L/2.L)/2.L, 3.L*std::sqrt(3.L/2.L)/2.L,
                                  std::sqrt(3.L/2.L),         -std::sqrt(3.L/2.L),    -3.L*std::sqrt(3.L/2.L)/2.L,
                                  -std::sqrt(3.L/2.L)/2.L;
            _innerRefCoeffs[2].engine().resize(6,0);
            _innerRefCoeffs[2] =                           -1.32957575266269475104L, 1.00928819006999249208L,
                                 1.11372017505655343851L,  -3.57158506454512406189L, 2.64084214163580360767L,
                                 0.614683388767643909165L;
            _innerRefCoeffs[3].engine().resize(6,0);
            _innerRefCoeffs[3] =                            3.35387951165860404196L, -0.141229992246845435862L,
                                -1.29954452099952139115L,  -0.599442481500795701774L, 1.08068229875787216245L,
                                 0.312647416473029867432L;
            _innerRefCoeffs[4].engine().resize(14,0);
            _innerRefCoeffs[4] =                           -0.000661624087242333507038L, 0.0191870985300276717041L,
                                -0.103410901756548045803L,  0.378884765460252640419L,   -1.07547200857749368479L,
                                 1.10969884455047509484L,   2.08619476135110666967L,     2.08619476135110666967L,
                                -0.795461007681955947706L,  0.000218300557375081050850L, 0.247205601265849432686L,
                                -0.146044023363594721918L,  0.0292055710131498807034L,  -0.00100708865562585795529L;
            _innerRefCoeffs[5].engine().resize(14,0);
            _innerRefCoeffs[5] =                           -0.00113465969390614975730L,  0.0329051311232783429617L,
                                -0.176288182431553100420L,  0.619104902646686873282L,   -1.70854321322203183450L,
                                 2.15119201578740835235L,   1.45657142326475405947L,    -3.40474081654023230539L,
                                 1.23029213497100169954L,  -0.000708279550556974276803L,-0.379632133141489874567L,
                                 0.224287719297057272008L, -0.0448526868857883735848L,   0.00154664437537201288224L;

            _innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[2] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[3] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[4] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[5] *= std::pow(2.L,-1.5L);

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  1;
            _innerOffsets[3] =  1;
            _innerOffsets[4] =  -7;
            _innerOffsets[5] =  -7;

            _innerH1SemiNorms = new long double[6];
            _innerH1SemiNorms[0] =  std::sqrt(13.9130434782608695652L);
            _innerH1SemiNorms[1] =  8.L;
            _innerH1SemiNorms[2] =  std::sqrt(268.675020609247651090L);
            _innerH1SemiNorms[3] =  std::sqrt(131.345247650197830689L);
            _innerH1SemiNorms[4] =  std::sqrt(81.2723006501163910509L);
            _innerH1SemiNorms[5] =  std::sqrt(263.708111955787835019L);

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _quadratic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(9,0);
            _rightSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(7,0);
            _rightRefCoeffs[0] =                              0.00142497525191845132977L, -0.0413242823056350885634L,
                                  0.221677722902601882504L,  -0.785751966578387326726L,    2.18220331826362106174L,
                                 -2.63492100553502467145L,   -2.39932206435886330782L;

            _rightRefCoeffs[0] *= std::pow(2.L,-1.5L);

            _rightOffsets = new long[1];
            _rightOffsets[0] =  1;

            _rightH1SemiNorms = new long double[1];
            _rightH1SemiNorms[0] =  std::sqrt(191.384235481065676847L);

            break;
            
        case 4:
            _addRefinementLevel = 2;
            _shiftFactor        = 8;
            // left part
            _numLeftParts = 5;
            _leftEvaluator = new Evaluator[5];
            _leftEvaluator[0] = _cubic_bspline_left_evaluator0;
            _leftEvaluator[1] = _cubic_bspline_inner_evaluator0;
            _leftEvaluator[2] = _cubic_bspline_inner_evaluator1;
            _leftEvaluator[3] = _cubic_bspline_inner_evaluator2;
            _leftEvaluator[4] = _cubic_bspline_inner_evaluator3;
            
            _leftSupport = new Support<T>[5];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            _leftSupport[4] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[5];
            _leftSingularSupport[0].engine().resize(5,0);
            _leftSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[1].engine().resize(3,0);
            _leftSingularSupport[1] = 0., 0.5, 1.;
            _leftSingularSupport[2].engine().resize(3,0);
            _leftSingularSupport[2] = 0., 0.5, 1.;
            _leftSingularSupport[3].engine().resize(5,0);
            _leftSingularSupport[3] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[4].engine().resize(5,0);
            _leftSingularSupport[4] = 0., 0.25, 0.5, 0.75, 1.;
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[5];
            _leftRefCoeffs[0].engine().resize(7,0);
            _leftRefCoeffs[0] =  7.53726229992987105824L,  -1.60433108517705699801L, -0.921289031788500302018L,
                                 1.68322742933961880394L,  -1.08863576448802049026L,  0.310369430432088613678L,
                                -0.0839025219232105910212L;
            _leftRefCoeffs[1].engine().resize(6,0);
            _leftRefCoeffs[1] =                                std::sqrt(35.L/13.L)/4.L, 3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L),       std::sqrt(35.L/13.L),     3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L)/4.L;
            _leftRefCoeffs[2].engine().resize(6,0);
            _leftRefCoeffs[2] =                                std::sqrt(35.L/3.L)/4.L,   std::sqrt(35.L/3.L)/2.L,
                                   std::sqrt(35.L/3.L)/2.L,    -std::sqrt(35.L/3.L)/2.L,  -std::sqrt(35.L/3.L)/2.L,
                                  -std::sqrt(35.L/3.L)/4.L;
            _leftRefCoeffs[3].engine().resize(6,0);
            _leftRefCoeffs[3] =                            -1.20968608645812591492L,   0.723769207408211462783L,
                                 2.78393856347267693448L,  -5.12943299583678560134L,   3.02163500540779053256L,
                                 0.583314706422121720912L;
            _leftRefCoeffs[4].engine().resize(6,0);
            _leftRefCoeffs[4] =                             3.23422528793847648562L,   0.295496503819005929826L,
                                 -2.64548042111713192471L,  0.336105140833428426303L,  0.872145478151140075300L,
                                  0.546845047769619692603L;

            _leftRefCoeffs[0] *= 0.5L;
            _leftRefCoeffs[1] *= 0.5L;
            _leftRefCoeffs[2] *= 0.5L;
            _leftRefCoeffs[3] *= 0.5L;
            _leftRefCoeffs[4] *= 0.5L;

            _leftOffsets = new long[5];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  1;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;
            _leftOffsets[4] =  1;

            _leftH1SemiNorms = new long double[5];
            _leftH1SemiNorms[0] =  std::sqrt(333.522181798281871759L);
            _leftH1SemiNorms[1] =  std::sqrt(12.9230769230769230769L);
            _leftH1SemiNorms[2] =  std::sqrt(56.L);
            _leftH1SemiNorms[3] =  std::sqrt(259.228164412374246759L);
            _leftH1SemiNorms[4] =  std::sqrt(123.272155615509401277L);


            // inner part
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_bspline_inner_evaluator0;
            _innerEvaluator[1] = _cubic_bspline_inner_evaluator1;
            _innerEvaluator[2] = _cubic_bspline_inner_evaluator2;
            _innerEvaluator[3] = _cubic_bspline_inner_evaluator3;
            _innerEvaluator[4] = _cubic_bspline_inner_evaluator4;
            _innerEvaluator[5] = _cubic_bspline_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(0.,1.);
            _innerSupport[3] = Support<T>(0.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 0.5, 1.;
            _innerSingularSupport[1].engine().resize(3,0);
            _innerSingularSupport[1] = 0., 0.5, 1.;
            _innerSingularSupport[2].engine().resize(5,0);
            _innerSingularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[3].engine().resize(5,0);
            _innerSingularSupport[3] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[4].engine().resize(9,0);
            _innerSingularSupport[4] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[5].engine().resize(9,0);
            _innerSingularSupport[5] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(6,0);
            _innerRefCoeffs[0] =                               std::sqrt(35.L/13.L)/4.L, 3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L),       std::sqrt(35.L/13.L),     3.L*std::sqrt(35.L/13.L)/4.L,
                                   std::sqrt(35.L/13.L)/4.L;
            _innerRefCoeffs[1].engine().resize(6,0);
            _innerRefCoeffs[1] =                                std::sqrt(35.L/3.L)/4.L,   std::sqrt(35.L/3.L)/2.L,
                                   std::sqrt(35.L/3.L)/2.L,    -std::sqrt(35.L/3.L)/2.L,  -std::sqrt(35.L/3.L)/2.L,
                                  -std::sqrt(35.L/3.L)/4.L;
            _innerRefCoeffs[2].engine().resize(6,0);
            _innerRefCoeffs[2] =                           -1.20968608645812591492L,   0.723769207408211462783L,
                                 2.78393856347267693448L,  -5.12943299583678560134L,   3.02163500540779053256L,
                                 0.583314706422121720912L;
            _innerRefCoeffs[3].engine().resize(6,0);
            _innerRefCoeffs[3] =                            3.23422528793847648562L,   0.295496503819005929826L,
                                 -2.64548042111713192471L,  0.336105140833428426303L,  0.872145478151140075300L,
                                  0.546845047769619692603L;
            _innerRefCoeffs[4].engine().resize(14,0);
            _innerRefCoeffs[4] =                            -0.00757430944378488466295L,  0.0304795243208619817595L,
                                 -0.172103404134940755815L,  0.560461227654867482281L,   -1.42578946715790175367L,
                                  0.970460974994275945818L,  2.18881605627051402000L,     2.18881605627051402000L,
                                 -0.619734802216407440776L, -0.418091821163284038735L,    0.735614181183121958188L,
                                 -0.475237518862456267198L,  0.135438721938371823978L,   -0.0366114708090683515524L;
            _innerRefCoeffs[5].engine().resize(14,0);
            _innerRefCoeffs[5] =                            -0.0183274146720283930429L,   0.0714158787346878889722L,
                                 -0.346389724181854150965L,  0.968551554773324101920L,   -2.38163361520215938629L,
                                  1.98214605626068896928L,   3.24419184855680683262L,    -4.62339059465912210039L,
                                  1.03257150206081408287L,   0.612554526218050167521L,   -1.11025696825143994205L,
                                  0.717898961851500034340L, -0.204656443589658444638L,    0.0553244321003903398448L;

            _innerRefCoeffs[0] *= 0.5L;
            _innerRefCoeffs[1] *= 0.5L;
            _innerRefCoeffs[2] *= 0.5L;
            _innerRefCoeffs[3] *= 0.5L;
            _innerRefCoeffs[4] *= 0.5L;
            _innerRefCoeffs[5] *= 0.5L;

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  1;
            _innerOffsets[3] =  1;
            _innerOffsets[4] = -7;
            _innerOffsets[5] = -7;

            _innerH1SemiNorms = new long double[6];
            _innerH1SemiNorms[0] =  std::sqrt(12.9230769230769230769L);
            _innerH1SemiNorms[1] =  std::sqrt(56.L);
            _innerH1SemiNorms[2] =  std::sqrt(259.228164412374246759L);
            _innerH1SemiNorms[3] =  std::sqrt(123.272155615509401277L);
            _innerH1SemiNorms[4] =  std::sqrt(79.7761595473175714506L);
            _innerH1SemiNorms[5] =  std::sqrt(258.536491525265347505L);

            // right parts
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(5,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(7,0);
            _rightRefCoeffs[0] =                              0.0227793146076243275116L,  -0.0890977886403828079147L,
                                  0.440560827664300141301L,  -1.25932003737469557193L,     3.11312728470113776594L,
                                 -2.51604382597377156721L,   -4.32606783787621990630L;

            _rightRefCoeffs[0] *= 0.5L;

            _rightOffsets = new long[1];
            _rightOffsets[0] =  1;

            _rightH1SemiNorms = new long double[1];
            _rightH1SemiNorms[0] =  std::sqrt(206.218396977329657205L);

            break;
            
        default: std::cerr << "BC for MRA<T,Orthogonal,Interval,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_MRA_TCC
