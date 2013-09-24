#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTMRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTMRA_TCC 1

#include <cassert>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
MRA<T,Orthogonal,Interval,MultiRefinement>::MRA(int _d, int j)
    : d(_d), j0((j<=0) ? 1 : j), _bc(2,0), _j(j0), phi(*this)
{
    assert(d>=1);
    assert(d==1 || j0>=1);
    if (d>1 && j0<1) {
        std::cerr << "MRA<T,Orthogonal,Interval,MultiRefinement>: "
                  << "d = " << d << " needs to be larger than 1. Stopping." << std::endl;
        exit(1);
    }

    /* L_2 orthonormal multiwavelet bases without Dirichlet boundary conditions are not
     * implemented yet. They require __different__ boundary adapted wavelets and scaling functions.
     */
    if (d==1) {
        // Same functions as in the standard mra. However, the computation of translation indices
        // is different in this class, so we need to change the arrangement of the evaluators.

        _shiftFactor        = 2;
        //left part
        _numLeftParts = 0;
        //_leftEvaluator = new Evaluator[0];
        //_leftSupport = new Support<T>[0];
        //_leftSingularSupport = new DenseVector<Array<T> >[0];
        //_leftRefCoeffs = new DenseVector<Array<long double> >[0];
        //_leftOffsets = new long[0];

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

        //right part
        _numRightParts = 0;
        //_rightEvaluator = new Evaluator[0];
        //_rightSupport = new Support<T>[0];
        //_rightSingularSupport = new DenseVector<Array<T> >[0];
        //_rightRefCoeffs = new DenseVector<Array<long double> >[0];
        //_rightOffsets   = new long[0];
    }
    else {
        this->enforceBoundaryCondition<DirichletBC>();
    }

    setLevel(_j);
}

template <typename T>
MRA<T,Orthogonal,Interval,MultiRefinement>::~MRA()
{
    if (_numLeftParts>0) {
        delete[] _leftEvaluator;
        delete[] _leftSupport;
        delete[] _leftSingularSupport;
        delete[] _leftRefCoeffs;
        delete[] _leftOffsets;
    }
    if (_numInnerParts>0) {
        delete[] _innerEvaluator;
        delete[] _innerSupport;
        delete[] _innerSingularSupport;
        delete[] _innerRefCoeffs;
        delete[] _innerOffsets;

    }
    if (_numRightParts>0) {
        delete[] _rightEvaluator;
        delete[] _rightSupport;
        delete[] _rightSingularSupport;
        delete[] _rightRefCoeffs;
        delete[] _rightOffsets;
    }
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Orthogonal,Interval,MultiRefinement>::cardI(int j) const
{
    assert(d==1 || j>=j0);
    // This is different to the corresponding multiscaling functions!!!
    // Case differentiation needed since piecewise linear cannot be smooth at the dyadic grid points
    if (d==1) return pow2i<int>(j);
    else if (d==2) return pow2i<int>(j)-1;
    else      return pow2i<int>(j)*(d-2);
}

template <typename T>
int
MRA<T,Orthogonal,Interval,MultiRefinement>::cardIL(int /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,MultiRefinement>::cardII(int j) const
{
    assert(d==1 || j>=j0);
    return cardI(j) - _numLeftParts - _numRightParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,MultiRefinement>::cardIR(int /*j*/) const
{
    return _numRightParts;
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,MultiRefinement>::rangeI(int j) const
{
    assert(d==1 || j>=j0);
    return Range<int>(0,cardI(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,MultiRefinement>::rangeIL(int /*j*/) const
{
    return Range<int>(0,cardIL() - 1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,MultiRefinement>::rangeII(int j) const
{
    assert(d==1 || j>=j0);
    return Range<int>(cardIL(), cardIL()+cardII(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,MultiRefinement>::rangeIR(int j) const
{
    assert(d==1 || j>=j0);
    return Range<int>(cardIL()+cardII(j),cardI(j)-1);
}

template <typename T>
int
MRA<T,Orthogonal,Interval,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,Interval,MultiRefinement>::setLevel(int j) const
{
    assert(d==1 || j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Orthogonal,Interval,MultiRefinement>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

    _numLeftParts = 0;
    _numRightParts = 0;
    switch (d) {
        case 2:
            _shiftFactor = 2;
            //left part
            _numLeftParts = 0;
            //_leftEvaluator = new Evaluator[0];
            //_leftSupport = new Support<T>[0];
            //_leftSingularSupport = new DenseVector<Array<T> >[0];
            //_leftRefCoeffs = new DenseVector<Array<long double> >[0];
            //_leftOffsets = new long[0];

            //inner part
            _numInnerParts = 1;
            _innerEvaluator = new Evaluator[1];
            _innerEvaluator[0] = _linear_refinement_inner_evaluator0;
            _innerSupport = new Support<T>[1];
            _innerSupport[0] = Support<T>(0.,2.);
            _innerSingularSupport = new DenseVector<Array<T> >[1];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 1., 2.;
            _innerRefCoeffs = new DenseVector<Array<long double> >[1];
            _innerRefCoeffs[0].engine().resize(3,0);
            _innerRefCoeffs[0] = 0.5L, 1.L, 0.5L;
            _innerRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _innerOffsets = new long[1];
            _innerOffsets[0] =  0;

            //right part
            _numRightParts = 0;
            //_rightEvaluator = new Evaluator[0];
            //_rightSupport = new Support<T>[0];
            //_rightSingularSupport = new DenseVector<Array<T> >[0];
            //_rightRefCoeffs = new DenseVector<Array<long double> >[0];
            //_rightOffsets = new long[0];

            break;

        case 3:
            _shiftFactor = 2;
            //left part
            _numLeftParts = 1;
            _leftEvaluator = new Evaluator[1];
            _leftEvaluator[0] = _quadratic_refinement_left_evaluator0;
            _leftSupport = new Support<T>[1];
            _leftSupport[0] = Support<T>(0.,2.);
            _leftSingularSupport = new DenseVector<Array<T> >[1];
            _leftSingularSupport[0].engine().resize(3,0);
            _leftSingularSupport[0] = 0., 1., 2.;

            _leftRefCoeffs = new DenseVector<Array<long double> >[1];
            _leftRefCoeffs[0].engine().resize(3,0);
            _leftRefCoeffs[0] = 0.5L, 0.75L, 0.25L;
            _leftRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _leftOffsets = new long[1];
            _leftOffsets[0] =  0;

            //inner part
            _numInnerParts = 1;
            _innerEvaluator = new Evaluator[1];
            _innerEvaluator[0] = _quadratic_refinement_inner_evaluator0;
            _innerSupport = new Support<T>[1];
            _innerSupport[0] = Support<T>(0.,3.);
            _innerSingularSupport = new DenseVector<Array<T> >[1];
            _innerSingularSupport[0].engine().resize(4,0);
            _innerSingularSupport[0] = 0., 1., 2., 3.;

            _innerRefCoeffs = new DenseVector<Array<long double> >[1];
            _innerRefCoeffs[0].engine().resize(4,0);
            _innerRefCoeffs[0] = 0.25L, 0.75L, 0.75L, 0.25L;
            _innerRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _innerOffsets = new long[1];
            _innerOffsets[0] =  1;

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _quadratic_refinement_right_evaluator0;
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,2.);
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(3,0);
            _rightSingularSupport[0] = 0., 1., 2.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(3,0);
            _rightRefCoeffs[0] = 0.25L, 0.75L, 0.5L;
            _rightRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _rightOffsets = new long[1];
            _rightOffsets[0] =  1;

            break;

        case 4:
            _shiftFactor = 4;
            //left part
            _numLeftParts = 1;
            _leftEvaluator = new Evaluator[1];
            _leftEvaluator[0] = _cubic_refinement_left_evaluator0;
            _leftSupport = new Support<T>[1];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSingularSupport = new DenseVector<Array<T> >[1];
            _leftSingularSupport[0].engine().resize(2,0);
            _leftSingularSupport[0] = 0., 1.;
            _leftRefCoeffs = new DenseVector<Array<long double> >[1];
            _leftRefCoeffs[0].engine().resize(3,0);
            _leftRefCoeffs[0] = 0.5L, 0.5L, 0.25L;
            _leftRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _leftOffsets = new long[1];
            _leftOffsets[0] =  0;

            //inner part
            _numInnerParts = 2;
            _innerEvaluator = new Evaluator[2];
            _innerEvaluator[0] = _cubic_refinement_inner_evaluator0;
            _innerEvaluator[1] = _cubic_refinement_inner_evaluator1;
            _innerSupport = new Support<T>[2];
            _innerSupport[0] = Support<T>(0.,2.);
            _innerSupport[1] = Support<T>(0.,2.);
            _innerSingularSupport = new DenseVector<Array<T> >[2];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 1., 2.;
            _innerSingularSupport[1].engine().resize(3,0);
            _innerSingularSupport[1] = 0., 1., 2.;
            _innerRefCoeffs = new DenseVector<Array<long double> >[2];
            _innerRefCoeffs[0].engine().resize(5,0);
            _innerRefCoeffs[0] = 0.25L, 0.625L, 0.75L, 0.25L, 0.125L;
            _innerRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _innerRefCoeffs[1].engine().resize(5,0);
            _innerRefCoeffs[1] = 0.125L, 0.25L, 0.75L, 0.625L, 0.25L;
            _innerRefCoeffs[1] *= std::pow(2.L,-0.5L);
            _innerOffsets = new long[2];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  2;

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_refinement_right_evaluator0;
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(1.,2.);
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(2,0);
            _rightSingularSupport[0] = 1., 2.;
            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(3,0);
            _rightRefCoeffs[0] = 0.25L, 0.5L, 0.5L;
            _rightRefCoeffs[0] *= std::pow(2.L,-0.5L);
            _rightOffsets = new long[1];
            _rightOffsets[0] =  5;

            break;


        default: std::cerr << "BSpline<T,Orthogonal,Interval,MultiRefinement> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_REFINEMENTMRA_TCC
