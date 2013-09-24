/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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
#include <numeric>
#include <cxxblas/cxxblas.h>
#include <extensions/flens/lapack.h>
#include <cxxblas/cxxblas.h>

namespace lawa {

// calculate values at integer nodes by solving an eigenvalue problem (EVP).
template <typename T>
void
_evalAtIntegersByEVP(const DenseVector<Array<T> > &a, 
                     DenseVector<Array<T> > &valuesAtIntegers)
{
    int l1 = a.firstIndex(),
        l2 = a.lastIndex();

    GeMatrix<FullStorage<T, cxxblas::ColMajor> > A(a.length(), a.length(), l1, l1);
    
    // fill matrix A of eigenvalue problem.
    for (int m=l1; m<=l2; ++m) {
        int from = std::max(l1, 2*m-l2),
              to = std::min(l2, 2*m-l1);
        for (int k=from; k<=to; ++k) {
            A(m,k) = a(2*m-k);
        }
    }
    
    // calculate eigenvalues of A.
    DenseVector<Array<T> > wr(A.numRows()), wi(A.numRows()); // eigenvalues
    GeMatrix<FullStorage<T,cxxblas::ColMajor> > VL, 
                                                VR(A.numRows(),A.numRows()); // left and right eigenvectors
    ev(false, true, A, wr, wi, VL, VR); //only calculate right eigenvectors.
    // TODO: as soon as iamax is in FLENS, take next line one level higher!!!!
    int pos = wr.firstIndex() + cblas_idamax(wr.length(), wr.engine().data(), 1);
    valuesAtIntegers = VR(_,pos); // choosing the corresponding eigenvector.
    
    // the elements of the eigenvector have to sum up to 1.
    T sum = std::accumulate(valuesAtIntegers.engine().data(),
                            valuesAtIntegers.engine().data()
                          + valuesAtIntegers.length(), 0.0);
    assert(sum);
    valuesAtIntegers /= sum;
    valuesAtIntegers.engine().changeIndexBase(l1);
}

template <typename T>
void
_evalAtIntegers(const BSpline<T,Primal,R,CDF> phi,
                DenseVector<Array<T> > &valuesAtIntegers)
{
    _evalAtIntegersByEVP(phi.a, valuesAtIntegers);
}

template <typename T>
void
_evalAtIntegers(const BSpline<T,Dual,R,CDF> phi_,
                DenseVector<Array<T> > &valuesAtIntegers)
{
    if (phi_.d==phi_.d_) {
        switch (phi_.d) {
            case 2: valuesAtIntegers.engine().resize(5,-2);
                    valuesAtIntegers = 0., -2., 5., -2., 0.;
                    
                    return;
            case 3: valuesAtIntegers.engine().resize(8,-3);
                    valuesAtIntegers = 0., -0.0075, -0.1025, 0.61,
                                       0.61, -0.1025, -0.0075, 0.;
                    return;
            case 4: valuesAtIntegers.engine().resize(11,-5);
                    valuesAtIntegers = 0.,
                                       0.00584608843537426023373448913389,
                                      -0.12627551020408175896925229153567,
                                      -0.28367346938775472864335824851878,
                                       2.01096938775509981311984120111447,
                                      -2.21373299319727800948953699844424,
                                       2.01096938775510247765510030149017,
                                      -0.28367346938775483966566071103443,
                                      -0.12627551020408162019137421339110,
                                       0.00584608843537434263309959803223,
                                      -0.0;
                    return;
           default: assert(0); // precalculated values not yet available!
                    return;
        }
    } else {
        _evalAtIntegersByEVP(phi_.a_, valuesAtIntegers);
    }
}

template <typename T>
void
subdivide(const BSpline<T,Primal,R,CDF> &phi, int j,
          DenseVector<Array<T> > &dyadicValues)
{
    DenseVector<Array<T> > valuesAtIntegers;
    _evalAtIntegers(phi, valuesAtIntegers);
    
    // set values at integer positions of result.
    int twoJ = pow2i<T>(j);
    int from = twoJ*phi.l1,
          to = twoJ*phi.l2;
    dyadicValues.engine().resize(to-from+1,from);
    dyadicValues(_(from,twoJ,to)) = valuesAtIntegers;
    
    // calculate values inbetween on dyadic grid.
    for (int l=1; l<=j; ++l) {
        for (int k=from+pow2i<T>(j-l); k<=to; k+=pow2i<T>(j-l+1)) {
            int mFrom = std::max(phi.l1, ((2*k-to)>>j)+1);
            int   mTo = std::min(phi.l2,  (2*k-from)>>j);
            for (int m=mFrom; m<=mTo; ++m) {
                dyadicValues(k) += phi.a(m)*dyadicValues(2*k-twoJ*m);
            }
        }
    }
}

template <typename T>
void
subdivide(const BSpline<T,Dual,R,CDF> &phi_, int j,
          DenseVector<Array<T> > &dyadicValues)
{
    DenseVector<Array<T> > valuesAtIntegers;
    _evalAtIntegers(phi_, valuesAtIntegers);
    
    // set values at integer positions of result.
    int twoJ = pow2i<T>(j);
    int from = twoJ*phi_.l1_,
          to = twoJ*phi_.l2_;
    dyadicValues.engine().resize(to-from+1,from);
    dyadicValues(_(from,twoJ,to)) = valuesAtIntegers;
    
    // calculate values inbetween on dyadic grid.
    for (int l=1; l<=j; ++l) {
        for (int k=from+pow2i<T>(j-l); k<=to; k+=pow2i<T>(j-l+1)) {
            int mFrom = std::max(phi_.l1_, ((2*k-to)>>j)+1);
            int   mTo = std::min(phi_.l2_,  (2*k-from)>>j);
            for (int m=mFrom; m<=mTo; ++m) {
                dyadicValues(k) += phi_.a_(m)*dyadicValues(2*k-twoJ*m);
            }
        }
    }
}

} // namespace lawa

