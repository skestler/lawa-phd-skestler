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
MRA<T,Primal,RPlus,SparseMulti>::MRA(int _d, int j)
    : d(_d), j0(j), _bc(1,0), _j(j0), phi(*this)
{
    assert(d>=2);

    setLevel(_j);
}

template <typename T>
MRA<T,Primal,RPlus,SparseMulti>::~MRA()
{
    delete[] _innerEvaluator;
    delete[] _innerSupport;
    delete[] _innerSingularSupport;
}

//--- cardinalities of left index sets. -------------------

template <typename T>
long
MRA<T,Primal,RPlus,SparseMulti>::cardIL(int /*j*/) const
{
    return _numLeftParts;
}

//--- ranges of left index sets. --------------------------

template <typename T>
Range<long>
MRA<T,Primal,RPlus,SparseMulti>::rangeIL(int /*j*/) const
{
    return Range<long>(0,0);
}

template <typename T>
int
MRA<T,Primal,RPlus,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,RPlus,SparseMulti>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Primal,RPlus,SparseMulti>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    switch (d) {
        case 4:
            // left B-splines
            _numLeftParts = 1;
            _leftEvaluator = new Evaluator[1];
            _leftEvaluator[0] = _cubic_sparsemulti_scaling_left_evaluator0;

            _leftSupport = new Support<T>[1];
            _leftSupport[0] = Support<T>(0,1);

            _leftSingularSupport = new DenseVector<Array<T> >[1];
            _leftSingularSupport[0] = linspace(0.,1.,2);

            _leftScalingFactors.engine().resize(1,0);
            _leftScalingFactors = 1.*std::sqrt(105.);

            // inner B-splines
            _numInnerParts = 2;
            _innerEvaluator = new Evaluator[2];
            _innerEvaluator[0] = _cubic_sparsemulti_scaling_inner_evaluator0;
            _innerEvaluator[1] = _cubic_sparsemulti_scaling_inner_evaluator1;

            _innerSupport = new Support<T>[2];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);

            _innerSingularSupport = new DenseVector<Array<T> >[2];
            _innerSingularSupport[0] = linspace(-1.,1.,3);
            _innerSingularSupport[1] = linspace(-1.,1.,3);

            _innerScalingFactors.engine().resize(2,0);
            _innerScalingFactors = std::sqrt(35./26.),std::sqrt(105./2.);


            break;

        default: std::cerr << "BSpline<T,Primal,RPlus,SparseMulti> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }



    _bc(0) = DirichletBC;

}

} // namespace lawa
