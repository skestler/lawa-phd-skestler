/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_wavelet_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Primal,Interval,SparseMulti>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _bc(2,0), _j(j0), _numSplines(4), psi(*this)
{
    if (d!=4) {
        std::cerr << "Basis<T,Primal,Interval,SparseMulti> not yet implemented for d=" << d << std::endl;
        exit(1);
    }
    if (d==4) {
        if (j0<1) {
            std::cerr << "Basis<T,Primal,Interval,SparseMulti> not yet implemented for j0=" << j0 << std::endl;
            exit(1);
        }
    }
    setLevel(_j);
}
    
template <typename T>
Basis<T,Primal,Interval,SparseMulti>::~Basis()
{
    delete[] _leftEvaluator;
    delete[] _innerEvaluator;
    delete[] _rightEvaluator;
    delete[] _leftSupport;
    delete[] _innerSupport;
    delete[] _rightSupport;
    delete[] _leftSingularSupport;
    delete[] _innerSingularSupport;
    delete[] _rightSingularSupport;
}

template <typename T>
int
Basis<T,Primal,Interval,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Interval,SparseMulti>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Primal,Interval,SparseMulti>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    _bc = 1,1;

    
    switch (d) {
        case 4:
            _numSplines = 4;        //required for lambdaTilde
            // left wavelets
            _numLeftParts = 3;
            _leftEvaluator = new Evaluator[3];
            _leftEvaluator[0] = _sparsemulti_cubic_wavelet_left_evaluator0;
            _leftEvaluator[1] = _sparsemulti_cubic_wavelet_left_evaluator1;
            _leftEvaluator[2] = _sparsemulti_cubic_wavelet_left_evaluator2;

            _leftSupport = new Support<T>[3];
            _leftSupport[0] = Support<T>(0.,2.);
            _leftSupport[1] = Support<T>(0.,2.);
            _leftSupport[2] = Support<T>(0.,2.);

            _leftSingularSupport = new DenseVector<Array<T> >[3];
            _leftSingularSupport[0] = linspace(0.0,2.0,5);
            _leftSingularSupport[1] = linspace(0.0,2.0,5);
            _leftSingularSupport[2] = linspace(0.0,2.0,5);

            _leftScalingFactors.engine().resize(3,0);
            _leftScalingFactors = 15./(2.*std::sqrt(13./7.)), 39./(8.*std::sqrt(11./7.)),
                                  60./std::sqrt(7841./1001.);

            // inner wavelets
            _numInnerParts = 4;
            _innerEvaluator = new Evaluator[4];
            _innerEvaluator[0] = _sparsemulti_cubic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _sparsemulti_cubic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _sparsemulti_cubic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _sparsemulti_cubic_wavelet_inner_evaluator3;

            _innerSupport = new Support<T>[4];
            _innerSupport[0] = Support<T> (0.,2.);
            _innerSupport[1] = Support<T> (0.,2.);
            _innerSupport[2] = Support<T>(-2.,2.);
            _innerSupport[3] = Support<T>(-2.,2.);

            _innerSingularSupport = new DenseVector<Array<T> >[4];
            _innerSingularSupport[0] = linspace(0.0,2.0,5);
            _innerSingularSupport[1] = linspace(0.0,2.0,5);
            _innerSingularSupport[2] = linspace(-2.0,2.0,9);
            _innerSingularSupport[3] = linspace(-2.0,2.0,9);

            _innerScalingFactors.engine().resize(4,0);
            _innerScalingFactors = 15./(2.*std::sqrt(13./7.)), 39./(8.*std::sqrt(11./7.)),
                                   30./std::sqrt(2467613./2002.), 30./std::sqrt(7841./2002.);

            // right wavelets
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _sparsemulti_cubic_wavelet_right_evaluator0;

            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(-2.0,0.0);

            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0] = linspace(-2.0,0.0,5);

            _rightScalingFactors.engine().resize(1,0);
            _rightScalingFactors = 60./std::sqrt(7841./1001.);

            break;

        default: std::cerr << "Wavelet<T,Primal,Interval,SparseMulti> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }


    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Primal,Interval,SparseMulti> &
Basis<T,Primal,Interval,SparseMulti>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

template <typename T>
Support<T>
Basis<T,Primal,Interval,SparseMulti>::max_support() const
{
    if (d==4) return Support<T>(-2.,2.);
}

template <typename T>
int
Basis<T,Primal,Interval,SparseMulti>::cardJ(int j) const
{
    assert(j>=j0);
    return 2*pow2i<int>(j);
}

template <typename T>
int
Basis<T,Primal,Interval,SparseMulti>::cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
int
Basis<T,Primal,Interval,SparseMulti>::cardJI(int j) const
{
    assert(j>=j0);
    return 2*pow2i<int>(j)-_numLeftParts-_numRightParts;
}

template <typename T>
int
Basis<T,Primal,Interval,SparseMulti>::cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Primal,Interval,SparseMulti>::rangeJ(int j) const
{
    assert(j>=j0);
    return Range<int>(1,cardJ(j));
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,SparseMulti>::rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<int>(1,cardJL(j));
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,SparseMulti>::rangeJI(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL(j)+1,cardJL(j)+cardJI(j));
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,SparseMulti>::rangeJR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL(j)+cardJI(j)+1,cardJ(j));
}




template <typename T>
long
Basis<T,Primal,Interval,SparseMulti>::long_cardJ(int j) const
{
    assert(j>=j0);
    return 2*pow2i<long>(j);
}

template <typename T>
long
Basis<T,Primal,Interval,SparseMulti>::long_cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
long
Basis<T,Primal,Interval,SparseMulti>::long_cardJI(int j) const
{
    assert(j>=j0);
    return 2*pow2i<long>(j)-_numLeftParts-_numRightParts;
}

template <typename T>
long
Basis<T,Primal,Interval,SparseMulti>::long_cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<long>
Basis<T,Primal,Interval,SparseMulti>::long_rangeJ(int j) const
{
    assert(j>=j0);
    return Range<long>(1,long_cardJ(j));
}

template <typename T>
const Range<long>
Basis<T,Primal,Interval,SparseMulti>::long_rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<long>(1,long_cardJL(j));
}

template <typename T>
const Range<long>
Basis<T,Primal,Interval,SparseMulti>::long_rangeJI(int j) const
{
    assert(j>=j0);
    return Range<long>(long_cardJL(j)+1,long_cardJL(j)+long_cardJI(j));
}

template <typename T>
const Range<long>
Basis<T,Primal,Interval,SparseMulti>::long_rangeJR(int j) const
{
    assert(j>=j0);
    return Range<long>(long_cardJL(j)+long_cardJI(j)+1,long_cardJ(j));
}

} // namespace lawa
