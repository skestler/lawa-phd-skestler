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

namespace lawa {

template <typename T>
Basis<T,Primal,RPlus,SparseMulti>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _bc(1,0), _j(j0), psi(*this)
{
    assert(d>=2);
    std::cerr << "Basis<T,Primal,RPlus,SparseMulti>::j0 = " << j0 << std::endl;
    setLevel(_j);
}
    
template <typename T>
Basis<T,Primal,RPlus,SparseMulti>::~Basis()
{
    delete[] _leftEvaluator;
    delete[] _innerEvaluator;
    delete[] _leftSupport;
    delete[] _innerSupport;
    delete[] _leftSingularSupport;
    delete[] _innerSingularSupport;
}

template <typename T>
int
Basis<T,Primal,RPlus,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,RPlus,SparseMulti>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Primal,RPlus,SparseMulti>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    _bc(0) = DirichletBC;

    
    switch (d) {
        case 4:
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

            break;

        default: std::cerr << "Wavelet<T,Primal,RPlus,SparseMulti> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }


    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Primal,RPlus,SparseMulti> &
Basis<T,Primal,RPlus,SparseMulti>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of left index sets (primal).
template <typename T>
long
Basis<T,Primal,RPlus,SparseMulti>::cardJL(int j) const
{
    assert(j>=j0);
    return _numLeftParts;
}

// ranges of whole, left index sets (primal).
template <typename T>
const Range<long>
Basis<T,Primal,RPlus,SparseMulti>::rangeJL(int j) const
{
    assert(j>=j0);
    return Range<long>(1,cardJL(j));
}

} // namespace lawa
