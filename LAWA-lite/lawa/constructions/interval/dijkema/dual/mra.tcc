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
#include <cmath>

#include <lawa/aux/aux.h>
#include <lawa/math/math.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <extensions/flens/lapack_flens.h>
#include <lawa/constructions/interval/primbs/primal/spline_helper.h>

namespace lawa {

template <typename T>
MRA<T,Dual,Interval,Dijkema>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), mu(d&1),
      l1((mu-d)/2), l2((mu+d)/2),
      l1_(l1-d_+1), l2_(l2+d_-1),
      min_j0(iceil<int>(log2(-2*l1_-1+mu))),
      j0((j==-1) ? min_j0 : j),
      phi_(*this),
      phi_R(_d,_d_),
      _bc(2,0), _j(j0)
{
    assert(d>1);
    assert((d+d_)%2==0);
    assert(_j>=min_j0);

    _calcM0_();
}

template <typename T>
MRA<T,Dual,Interval,Dijkema>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Dual,Interval,Dijkema>::cardI_(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
int
MRA<T,Dual,Interval,Dijkema>::cardI_L(int /*j*/) const
{
    return d + d_ - 2 -_bc(0);
}

template <typename T>
int
MRA<T,Dual,Interval,Dijkema>::cardI_I(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) - 2*(d + d_ - 2);
}

template <typename T>
int
MRA<T,Dual,Interval,Dijkema>::cardI_R(int /*j*/) const
{
    return d + d_ - 2 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_(int j) const
{
    assert(j>=min_j0);
    return Range<int>(1 + _bc(0), pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_L(int /*j*/) const
{
    return Range<int>(1 + _bc(0), d + d_ - 2);
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_I(int j) const
{
    assert(j>=min_j0);
    return Range<int>(d + d_ - 1, pow2i<T>(j) - d_ + 1);
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_R(int j) const
{
    assert(j>=min_j0);
    return Range<int>(pow2i<T>(j) - d_ + 2, pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
int
MRA<T,Dual,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M0_.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Dual,Interval,Dijkema>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

    _calcM0_();
}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::_calcM0_()
{
    typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;
    typedef DenseVector<Array<T> > DenseVector;
    
    const int cons_j = ((d==2) && ((d_==2) || (d_==4))) ? min_j0+1 : min_j0;

    int extra = std::max(2-d+_bc(0),0);
    FullColMatrix A(_(-l2_+1,-l1_-1+extra),_(-l2_+1,-l1_-1+extra));
    for (int k=A.firstRow(); k<=A.lastRow(); ++k) {
        for (int m=A.firstCol(); m<=A.lastCol(); ++m) {
            if (l1_<=k-2*m && k-2*m<=l2_) {
                A(k,m) = phi_R.a_(k-2*m);
            }
        }
    }

    FullColMatrix VL, VR(A.numRows(),A.numRows());
    DenseVector r(VR.numRows()), i(VR.numRows());
    ev(false, true, A, r, i, VL, VR);

    DenseVector r_sorted = r;
    std::sort(r_sorted.engine().data(),
              r_sorted.engine().data()+r_sorted.length(),
              descending_by_magnitude<T>());

    DenseVector positions(r.length());
    for (int i=r.firstIndex(); i<=r.lastIndex(); ++i) {
        int j=0;
        while (r(i)!=r_sorted(++j)) { ; }
        positions(i) = j;
    }

    DenseVector indices;
    if ((_bc(0)==1) && (_bc(1)==1)) {
        if ((d==3) && (d_==3)) {
            indices.engine().resize(3);
            indices = 1, 2, 1;
        } else if ((d==3) && (d_==5)) {
            indices.engine().resize(5);
            indices = 1, 3, 6, 8, 10;
        } else if ((d==3) && (d_==7)) {
            indices.engine().resize(7);
            indices = 1, 2, 5, 8, 10, 12, 14;
        } else if ((d==3) && (d_==9)) {
            indices.engine().resize(9);
            indices = 1, 2, 5, 7, 10, 12, 14, 15, 18;
        } else if ((d==4) && (d_==6)) {
            indices.engine().resize(7);
            indices = 1, 4, 6, 7, 9, 11, 13;
        } else if ((d==4) && (d_==8)) {
            indices.engine().resize(9);
            indices = 1, 2, 4, 6, 9, 11, 13, 14, 17;
        } else if ((d==4) && (d_==10)) {
            indices.engine().resize(11);
            indices = 1, 2, 5, 8, 10, 12, 13, 16, 17, 18, 21;
        } else if ((d==5) && (d_==9)) {
            indices.engine().resize(11);
            indices = 1, 2, 4, 7, 9, 12, 13, 15, 16, 17, 20;
        } else if ((d==2) && (d_==2)) {
            indices.engine().resize(2);
            //--------------------------
            //indices = 1, 1;
            indices = 2, 4;
            positions = 1, 2, 3, 4;
            //--------------------------
        } else if ((d==2) && (d_==4)) {
            indices.engine().resize(4);
            indices = 1, 1, 4, 6;
        } else if ((d==2) && (d_==6)) {
            indices.engine().resize(6);
            indices = 1, 6, 8, 12, 9, 11;
        } else if ((d==2) && (d_==8)) {
            indices.engine().resize(8);
            indices = 1, 5, 4, 16, 11, 12, 9, 13;
        } else { // for d=2 we need all, for others not necessarily optimal!
            indices.engine().resize(d+d_-2);
            for (int i=indices.firstIndex(); i<=indices.lastIndex(); ++i) {
                indices(i) = i;
            }
        }
    } else {
        assert(_bc(0)==0);
        assert(_bc(1)==0);
        if ((d==3) && (d_==3)) {
            indices.engine().resize(4);
            indices = 2,4,5,6;
        } else if ((d==3) && (d_==5)) {
            indices.engine().resize(6);
            indices = 1, 3, 6, 8, 10, 4;
        } else if ((d==3) && (d_==7)) {
            indices.engine().resize(8);
            indices = 1, 2, 5, 7, 8, 10, 12, 14;
        } else if ((d==3) && (d_==9)) {
            indices.engine().resize(10);
            indices = 1, 2, 3, 5, 7, 10, 12, 14, 15, 18;
        } else if ((d==4) && (d_==6)) {
            indices.engine().resize(8);
            indices = 1, 4, 7, 9, 11, 13, 3, 8;
        } else if ((d==4) && (d_==8)) {
            indices.engine().resize(10);
            indices = 1, 3, 4, 5, 6, 9, 11, 13, 14, 17;
        } else if ((d==4) && (d_==10)) {
            indices.engine().resize(12);
            indices = 1, 2, 5, 7, 8, 9, 10, 13, 16, 17, 18, 21;
        } else if ((d==5) && (d_==9)) {
            indices.engine().resize(12);
            indices = 1, 3, 4, 5, 7, 9, 11, 12, 15, 16, 17, 20;
            //--------------------------
        } else if ((d==2) && (d_== 2)) {
            indices.engine().resize(2);
            indices = 2, 1;
            positions = 1, 2, 3;
            //--------------------------
        } else if ((d==2) && (d_==4)) {
            indices.engine().resize(4);
            indices = 1, 4, 5, 2;
        } else if ((d==2) && (d_==2)) {
            indices.engine().resize(2);
            indices = 2,1;
            positions = 1,2,3;
        } else { // for d=2 we need all, for others not necessarily optimal!
            indices.engine().resize(d+d_-2);
            for (int i=indices.firstIndex(); i<=indices.lastIndex(); ++i) {
                indices(i) = i;
            }
        }
    }

    FullColMatrix C_L(d+2*d_-3+extra, -l1_-1+l2+extra-_bc(0)-_bc(1));
    for (int i=1+_bc(0); i<=indices.length(); ++i) {
        C_L(_,i-_bc(0)) = VR(_,positions(indices(i)));
    }

    FullColMatrix C_R;
    arrow(C_L,C_R);

    FullColMatrix R_init(pow2i<T>(cons_j)+l2_-l1_-1,
                         pow2i<T>(cons_j)+l2-l1-1-_bc(0)-_bc(1),
                         -l2_+1, -l2+1+_bc(0));
    for (int c=R_init.firstCol(); c<=R_init.lastCol(); ++c) {
        R_init(c,c) = 1.;
    }
    C_L.engine().changeIndexBase(-l2_+1, -l2+1+_bc(0));
    R_init(C_L) = C_L;
    C_R.engine().changeIndexBase(R_init.lastRow()-C_R.numRows()+1,
                                 R_init.lastCol()-C_R.numCols()+1);
    R_init(C_R) = C_R;
    BSpline<T,Primal,R,CDF> phi(d);
    int kmin = -l2+1, mmin = -l2_+1, pmin = -l2+1, qmin = -l2_+1,
        kmax = -l1-1, mmax = -l1_-1, pmax = -l1-1, qmax = -l1_-1;
    FullColMatrix C(_(pmin,pmax), _(qmin,qmax));
    int ZLength = C.numRows()*C.numCols();
    FullColMatrix Z(_(pmin*(qmax-qmin+1)+qmin,
                      pmin*(qmax-qmin+1)+qmin+ZLength-1),
                    _(pmin*(qmax-qmin+1)+qmin, 
                      pmin*(qmax-qmin+1)+qmin+ZLength-1));
    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            for (int k=kmin; k<=kmax; ++k) {
                for (int m=mmin; m<=mmax; ++m) {
                    if (l1<=k-2*p && k-2*p<=l2 &&
                        l1_<=m-2*q && m-2*q<=l2_) {
                            Z(p*(qmax-qmin+1)+q, k*(mmax-mmin+1)+m) =
                                           0.5 * phi.a(k-2*p) * phi_R.a_(m-2*q);
                    }
                }
            }
        }
    }
    for (int i=Z.firstRow(); i<=Z.lastRow(); ++i) {
        Z(i,i) -= 1.;
    }
    Z *= -1;

    DenseVector f(Z.numCols(), Z.firstCol());
    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            T sum = 0.;
            for (int k=std::max(-l2+1,l1+2*p); k<=l2+2*p; ++k) {
                for (int m=std::max(std::max(-l2_+1,l1_+2*q),l1-l2_+1+k);
                     m<=std::min(l2_+2*q, l2-l1_-1+k); ++m) {
                    if (!((kmin<=k && k<=kmax) && (mmin<=m && m<=mmax))) {
                        sum += 0.5 * phi.a(k-2*p) * phi_R.a_(m-2*q) * (k==m);
                    }
                }
            }
            f(p*(qmax-qmin+1)+q) = sum;
        }
    }
    flens::DenseVector<Array<int> > pivots(f.length(),f.firstIndex());
    sv(Z,pivots,f);
    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            C(p,q) = f(p*(qmax-qmin+1)+q);
        }
    }
    FullColMatrix ExtMass(pow2i<T>(cons_j)+d-1,
                          pow2i<T>(cons_j)+l2_-l1_-1,
                          -l2+1, -l2_+1);
    for (int r=ExtMass.firstRow(); r<=ExtMass.lastRow(); ++r) {
        ExtMass(r,r) = 1.;
    }
    ExtMass(C) = C;
    FullColMatrix CR;
    arrow(C,CR);
    CR.engine().changeIndexBase(ExtMass.lastRow()-C.numRows()+1,
                                ExtMass.lastCol()-C.numCols()+1);
    ExtMass(CR) = CR;

    //--- duplicate code from primal Primbs-MRA method _calc_M0 ----------------
    FullColMatrix R(_(1-l2, pow2i<T>(cons_j)-l1-1),
                    _(1-l2+_bc(0), pow2i<T>(cons_j)-l1-1-_bc(1)));
    int rr = R.firstRow()+_bc(0);
    for (int c=R.firstCol(); c<=R.lastCol(); ++c, ++rr) {
        R(rr,c) = 1.;
    }
    // replace this part as soon as well understood ;-)
    DenseVector knots = linspace(-d+1., d-1., 2*d-1);
    knots.engine().changeIndexBase(1);
    FullColMatrix Transformation(knots.length()-d, knots.length()-d);
    Transformation.diag(0) = 1.;
    for (int i=1; i<d; ++i) {
        FullColMatrix Tmp = insertKnot(d-1,knots,(T)0.), Tmp2;
        blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,
                 1.,Tmp,Transformation,0.,Tmp2);
        Transformation = Tmp2;
    }

    FullColMatrix InvTrans = Transformation(_(Transformation.firstRow()+(d-1),
                                              Transformation.lastRow()),
                                            _(Transformation.firstCol(),
                                              Transformation.lastCol()) );

    //--- inverse(InvTrans)
    FullColMatrix TransTmp = InvTrans, Trans, TransTmp2;
    Trans = inv(InvTrans);
    
    Trans.engine().changeIndexBase(R.firstRow(),R.firstCol());
    if (d>2) {
        R(_(Trans.firstRow(),Trans.lastRow()),
          _(Trans.firstCol(),Trans.lastCol()-_bc(0)))
              = Trans( _ , _(Trans.firstCol()+_bc(0),Trans.lastCol()));

        arrow(Trans,TransTmp);
        TransTmp.engine().changeIndexBase(R.lastRow()-TransTmp.numRows()+1,
                                          R.lastCol()-TransTmp.numCols()+1);
        R(_(TransTmp.firstRow(),TransTmp.lastRow()),
          _(TransTmp.firstCol()+_bc(1),TransTmp.lastCol()))
              = TransTmp( _ , _(TransTmp.firstCol(),TransTmp.lastCol()-_bc(1)));
    }

    FullColMatrix ExtM0(_(-l2+1,pow2i<T>(cons_j+1)-l1-1),
                        _(-l2+1,pow2i<T>(cons_j)-l1-1));
    for (int q=ExtM0.firstCol(); q<=ExtM0.lastCol(); ++q) {
        for (int p=std::max(l1+2*q,ExtM0.firstRow());
             p<=std::min(l2+2*q,ExtM0.lastRow()); ++p) {
            ExtM0(p,q) = phi.a(p-2*q);
        }
    }

    FullColMatrix RjPlus1(_(-l2+1,pow2i<T>(cons_j+1)-l1-1),
                          _(-l2+1+_bc(0),pow2i<T>(cons_j+1)-l1-1-_bc(1)));
    rr = RjPlus1.firstRow()+_bc(0);
    for (int c=RjPlus1.firstCol(); c<=RjPlus1.lastCol(); ++c, ++rr) {
        RjPlus1(rr,c) = 1.;
    }
    if (d>2) {
        Trans.engine().changeIndexBase(RjPlus1.firstRow(),RjPlus1.firstCol());
        RjPlus1(_(Trans.firstRow(),Trans.lastRow()),
                _(Trans.firstCol(),Trans.lastCol()-_bc(0)))
            = Trans( _ , _(Trans.firstCol()+_bc(0), Trans.lastCol()));
        TransTmp.engine().changeIndexBase(
                                         RjPlus1.lastRow()-TransTmp.numRows()+1,
                                         RjPlus1.lastCol()-TransTmp.numCols()+1
                                         );
        RjPlus1(_(TransTmp.firstRow(),TransTmp.lastRow()),
                _(TransTmp.firstCol()+_bc(1),TransTmp.lastCol()))
            = TransTmp( _ , _(TransTmp.firstCol(), TransTmp.lastCol()-_bc(1)));
    }
    //--------------------------------------------------------------------------
    FullColMatrix Mass, MassTmp;
    blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
             1., ExtMass, R_init, 0., MassTmp);

    blas::mm(cxxblas::Trans, cxxblas::NoTrans,
             1., R, MassTmp, 0., Mass);

    FullColMatrix InvMass, InvMassTmp = Mass;
    InvMass = inv(Mass);

    FullColMatrix R_(R_init.numRows(), R_init.numCols(),
                     R_init.firstRow(), R_init.firstCol());
    blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1., R_init, InvMass, 0., R_);

    R_Left = R_(_(C_L.firstRow(),C_L.lastRow()),
                _(C_L.firstCol(),C_L.lastCol()+_bc(0)));
    R_Right = R_(_(C_R.firstRow(),C_R.lastRow()),
                 _(C_R.firstCol()-_bc(1),C_R.lastCol()));

    FullColMatrix ExtM0_(_(-l2_+1,pow2i<T>(cons_j+1)-l1_-1),
                         _(-l2_+1,pow2i<T>(cons_j)-l1_-1));
    for (int q=ExtM0_.firstCol(); q<=ExtM0_.lastCol(); ++q) {
        for (int p=std::max(l1_+2*q, ExtM0_.firstRow());
             p<=std::min(l2_+2*q, ExtM0_.lastRow()); ++p) {
            ExtM0_(p,q) = phi_R.a_(p-2*q);
        }
    }

    FullColMatrix ExtMassjPlus1(pow2i<T>(cons_j+1)+d-1, 
                                pow2i<T>(cons_j+1)+l2_-l1_-1,
                                -l2+1, -l2_+1);
    for (int r=ExtMassjPlus1.firstRow(); r<=ExtMassjPlus1.lastRow(); ++r) {
        ExtMassjPlus1(r,r) = 1.;
    }
    ExtMassjPlus1(C) = C;
    CR.engine().changeIndexBase(ExtMassjPlus1.lastRow()-C.numRows()+1,
                                ExtMassjPlus1.lastCol()-C.numCols()+1);
    ExtMassjPlus1(CR) = CR;
    FullColMatrix Mj0_, M0_Tmp;
    blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,RjPlus1,ExtMassjPlus1,0.,Mj0_);
    blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1., Mj0_, ExtM0_, 0., M0_Tmp);
    blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1., M0_Tmp, R_, 0., Mj0_);
    Mj0_.engine().changeIndexBase(1+_bc(0),1+_bc(1));
    blas::scal(Const<T>::R_SQRT2, Mj0_);
    R_Left.engine().changeIndexBase(1,1+_bc(0));
    R_Right.engine().changeIndexBase(1,1);
    M0_ = RefinementMatrix<T,Interval,Dijkema>(d+d_-2-_bc(0), d+d_-2-_bc(1),
                                               Mj0_, min_j0, cons_j);
    setLevel(_j);
}

} // namespace lawa

