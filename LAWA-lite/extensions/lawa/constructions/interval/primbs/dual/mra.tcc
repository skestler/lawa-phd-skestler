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
#include <list>

#include <lawa/aux/arrow.h>
#include <lawa/math/math.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <extensions/flens/lapack_flens.h>

namespace lawa {

template <typename T>
MRA<T,Dual,Interval,Primbs>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), mu(d&1),
      min_j0(iceil<int>(log(d+2*d_-3)/log(2))+1),
      j0((j==-1) ? min_j0 : j), phi_R(d,d_), phi_(),
      _bc(2,0), _j(j0)
{
    assert(d>1);
    assert(_j>=min_j0);

    _calcM0_();
}

template <typename T>
MRA<T,Dual,Interval,Primbs>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Dual,Interval,Primbs>::cardI_(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
int
MRA<T,Dual,Interval,Primbs>::cardI_L(int /*j*/) const
{
    return d + d_ - 2 -_bc(0);
}

template <typename T>
int
MRA<T,Dual,Interval,Primbs>::cardI_I(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - 2*(d+d_-2);
}

template <typename T>
int
MRA<T,Dual,Interval,Primbs>::cardI_R(int /*j*/) const
{
    return d + d_ - 2 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Dual,Interval,Primbs>::rangeI_(int j) const
{
    assert(j>=min_j0);
    return Range<int>(1 + _bc(0), pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Primbs>::rangeI_L(int /*j*/) const
{
    return Range<int>(1 + _bc(0), d + d_ - 2);
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Primbs>::rangeI_I(int j) const
{
    assert(j>=min_j0);
    return Range<int>(d+d_-1, pow2i<T>(j)+d-1-(d+d_-2));
}

template <typename T>
Range<int>
MRA<T,Dual,Interval,Primbs>::rangeI_R(int j) const
{
    assert(j>=min_j0);
    return Range<int>(pow2i<T>(j)+d-1-(d+d_-2)+1, pow2i<T>(j)+d-1-_bc(1));
}

template <typename T>
int
MRA<T,Dual,Interval,Primbs>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::setLevel(int j) const
{
    if (j!=_j) {
        assert(j>=min_j0);
        _j = j;
        M0_.setLevel(_j);
    }
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Dual,Interval,Primbs>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

    _twoScaleDual_2(d, d_);

    GeMatrix<FullStorage<T,ColMajor> > A(d+d_-2,d+d_-2);
    for (int n=1; n<=d+d_-2; ++n) {
        for (int m=1; m<=d+d_-2; ++m) {
            A(n,m) = MDD(m,n);
        }
    }

    for (int n=1; n<=d+d_-2; ++n) {
        A(n,n) = A(n,n)-1;
    }
    DenseVector<Array<T> > O(d+d_-2), x, tmp;
    A(d+d_-2,d+d_-2) = 1.;
    O(d+d_-2) = 1.;

    x = inv(A)*O;
    _orthonormalize(MP,MDD);
    tmp = x;
    x = transpose(inv(TT))*tmp;
    x.engine().changeIndexBase(1);
    GeMatrix<FullStorage<T,ColMajor> > TR(d+d_-2,d+d_-2);
    for (int n=1; n<=d+d_-2; ++n) {
        TR(n,n) = 1;
        if (n>1) {
            TR(n,1) = -x(n) / x(1);
        }
    }

    GeMatrix<FullStorage<T,ColMajor> > TRD(3*d_+2*d-5,3*d_+2*d-5), MDDD_H, Tmp, InvTR;
    TRD.diag(0) = 1.;

//    TRD(_(1,d+d_-2),_(1,d+d_-2)) = transpose(inv(TR));
    InvTR = inv(TR);
    for (int r=1; r<=d+d_-2; ++r) {
        for (int c=1; c<=d+d_-2; ++c) {
            TRD(r,c) = InvTR(c,r);
        }
    }

    blas::mm(NoTrans,NoTrans,1.,TRD,MDDD,0.,Tmp);
    blas::mm(NoTrans,Trans,1.,Tmp,TR,0.,MDDD_H);

    T factor = Const<T>::R_SQRT2;

    GeMatrix<FullStorage<T,ColMajor> > Mj0_(pow2i<T>(min_j0+1)+d-3,
                                            pow2i<T>(min_j0)+d-3, 
                                            2,2),
                                       Mj0_Right, Bugfix;
    Bugfix = MDDD_H(_(2,MDDD.numRows()), _(2,MDDD.numCols()));
    blas::scal(factor,Bugfix);
    Mj0_(_(2,MDDD.numRows()), _(2,MDDD.numCols())) = Bugfix;
    int row = d+d_-1;
    for (int c=MDDD.lastCol()+1; c<=Mj0_.lastCol()-MDDD.numCols()+1; ++c,row+=2) {
        Mj0_(_(row,row+phi_R.a_.length()-1),c) = factor*phi_R.a_;
    }
    arrow(Bugfix,Mj0_Right);
    Mj0_(_(Mj0_.lastRow()-Mj0_Right.numRows()+1,Mj0_.lastRow()),
         _(Mj0_.lastCol()-Mj0_Right.numCols()+1,Mj0_.lastCol())) = Mj0_Right;
    M0_ = RefinementMatrix<T,Interval,Primbs>(Mj0_Right.numCols(),Mj0_Right.numCols(), 
                                              Mj0_, min_j0, min_j0);
    M0_.setLevel(_j);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_orthonormalize(const GeMatrix<FullStorage<T,ColMajor> > &L3,
                                             const GeMatrix<FullStorage<T,ColMajor> > &LL3)
{
    int n = L3.numRows();
    int m = L3.numCols();

    GeMatrix<FullStorage<T,ColMajor> > L2(n-m,m), LL2(n-m,m), L22, C(m,m);

    for (int j=1; j<=n-m; ++j) {
        for (int k=1; k<=m; ++k) {
            L2(j,k) = L3(j+m,k);
            LL2(j,k) = LL3(j+m,k);
        }
    }

    // L22 = transpose(L2);
    L22.engine().resize(L2.cols(),L2.rows());
    for (int r=L2.firstRow(); r<=L2.lastRow(); ++r) {
        for (int c=L2.firstCol(); c<=L2.lastCol(); ++c) {
            L22(c,r) = L2(r,c);
        }
    }
    
    blas::mm(NoTrans,NoTrans,0.5,L22,LL2,0.,C);
    DenseVector<Array<T> > c(powii(m,2));

    for (int j=1; j<=m; ++j) {
        for (int k=1; k<=m; ++k) {
            c((j-1)*m+k) = C(j,k);
        }
    }
    GeMatrix<FullStorage<T,ColMajor> > B(powii(m,2),powii(m,2));

    for (int j=1; j<=m; ++j) {
        for (int k=1; k<=m; ++k) {
            for (int l=1; l<=m; ++l) {
                for (int s=1; s<=m; ++s) {
                    B((j-1)*m+k,(l-1)*m+s) = -0.5*L3(l,j)*LL3(s,k);
                }
            }
        }
    }

    for (int j=1; j<=powii(m,2); ++j) {
        B(j,j) += 1;
    }

    DenseVector<Array<T> > t(powii(m,2)), cc, dd;
    bool singular;
    B.engine().changeIndexBase(0,0);
    c.engine().changeIndexBase(0);
    qrf(B,cc,dd,singular);
    qrsolv(B,cc,dd,c);
    t = c;
    t.engine().changeIndexBase(1);
    TT.engine().resize(m,m);

    for (int j=1; j<=m; ++j) {
        for (int k=1; k<=m; ++k) {
            TT(j,k) = t((j-1)*m+k);
        }
    }
}

//------------------------------------------------------------------------------

template <typename T>
T
_dividedDifferences(DenseVector<Array<T> > &a, T x)
{
    int m = a.length();
    DenseVector<Array<T> > b(m);

    for (int j=1; j<=m; ++j) {
        b(j) = std::pow(std::max(a(j)-x,T(0)),m-2);
    }

    for (int j=1; j<=m-1; ++j) {
        for (int k=1; k<=m-j; ++k) {
            if (a(k)==a(k+j)) {
                b(k) = 0;
            } else {
                b(k) = (b(k+1)-b(k)) / (a(k+j)-a(k));
            }
        }
    }
    return b(1);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_Umrechnung1(int d)
{
    GeMatrix<FullStorage<T,ColMajor> > B(d-1,d-1), C(d-1,d-1);
    int m = d+1;
    DenseVector<Array<T> > b(m), c(m);
    for (int j=1; j<=m; ++j) {
        c(j) = j-m;
    }
    for (int l=1; l<=d-1; ++l) {
        for (int j=1; j<=m-1; ++j) {
            b(j) = b(j+1); c(j) = c(j)+1;
        }
        ++b(m); ++c(m);
        for (int i=1; i<=d-1; ++i) {
            B(i,l) = l*_dividedDifferences(b,T(i-1));
            C(i,l) = d*_dividedDifferences(c,T(i-1));
        }
    }
    blas::mm(NoTrans,NoTrans,1.,inv(C),B,0.,W);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_Umrechnung2(int d)
{
    DenseVector<Array<T> > a(2*d+1);
    GeMatrix<FullStorage<T,ColMajor> > Tmp, InvW, B(2*(d-1),2*(d-1));
    Mj.engine().resize(2*(d-1),d-1);

    for (int j=1; j<=d+1; ++j) {
        a(2*d+2-j) = pow2i<T>(1-d)*binomial(d,j-1);
    }

    for (int j=1; j<=d-1; ++j) {
        for (int k=1; k<=2*j; ++k) {
            Mj(k,j) = a(2*d+1-2*j+k);
        }
    }

    _Umrechnung1(d);
    InvW = inv(W);

    for (int j=1; j<=d-1; ++j) {
        for (int k=1; k<=d-1; ++k) {
            B(j,k) = InvW(j,k);
        }
    }

    for (int j=d; j<=2*(d-1); ++j) {
        B(j,j) = 1;
    }
    blas::mm(NoTrans,NoTrans,1.,B,Mj,0.,Tmp);
    blas::mm(NoTrans,NoTrans,1.,Tmp,W,0.,Mp);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScalePrimal(int d, int d_)
{
    MP.engine().resize(3*d_+2*d-5, d+d_-2);
    _Umrechnung2(d);

    for (int j=1; j<=2*d-2; ++j) {
        for (int k=1; k<=d-1; ++k) {
            MP(j,k) = Mp(j,k);
        }
    }

    DenseVector<Array<T> > b = BSpline<T,Primal,R,CDF>(d).a;
    b.engine().changeIndexBase(1);
    for (int k=d; k<=d+d_-2; ++k) {
        for (int j=2*k-d; j<=2*k; ++j) {
            MP(j,k) = b(j-2*k+d+1);
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScaleDual_1(int d, int d_)
{
    M1.engine().resize(3*d_+2*d-5,d+d_-2) || M1.engine().fill();
    MD.engine().resize(3*d_+2*d-5,d+d_-2) || MD.engine().fill();
    DenseVector<Array<T> > b = BSpline<T,Primal,R,CDF>(d).a;
    DenseVector<Array<T> > a = BSpline<T,Dual,R,CDF>(d,d_).a_;
    b.engine().changeIndexBase(1); a.engine().changeIndexBase(1);

    for (int l=d/2+d_-1; l<=d/2+3*d_+d-5; ++l) {
        int gu = std::max(0,iceil<int>((d/2+d_-3-l)/2.));
        int go = std::min(2*d_+d-4, (d+d/2+3*d_-5-l)/2);
        for (int n=gu; n<=go; ++n) {
            M1(l+iceil<int>(d/2.),d+d_-2) += a(l+2*n-d/2-d_+4);
        }
        for (int k=1; k<=d_-1; ++k) {
            for (int n=gu; n<=go; ++n) {
                M1(l+iceil<int>(d/2.),d+d_-2-k) += powii(n,k)*a(l+2*n-d/2-d_+4);
            }
        }
    }

    DenseVector<Array<T> > m(d_);
    for (int n=-1; n<=d+2*d_-3; ++n) {
        m(1) += a(n+2);
    }
    for (int j=2; j<=d_; ++j) {
        for (int n=-1; n<=d+2*d_-3; ++n) {
            m(j) += a(n+2)*powii(n,j-1);
        }
    }
    m *= .5;

    for (int k=0; k<=d_-1; ++k) {
        for (int r=0; r<=k; ++r) {
            M1(d+d_-2-r,d+d_-2-k) = pow2i<T>(-k)*binomial(k,r)*m(k-r+1);
        }
    }

    GeMatrix<FullStorage<T,ColMajor> > D3(d_,d_);
    D3(1,1) = 1;
    for (int k=2; k<=d_; ++k) {
        for (int i=1; i<=k; ++i) {
            for (int s=0; s<=i-1; ++s) {
                D3(k,i) += powii(-1,(i-1)-s)*binomial(i-1,s)*powii(s,k-1);
            }
        }
    }

    GeMatrix<FullStorage<T,ColMajor> > J(d_,d_);
    for (int j=1; j<=d_; ++j) {
        J(d_+1-j,j) = 1;
    }

    GeMatrix<FullStorage<T,ColMajor> > D4, Tmp;
    blas::mm(NoTrans,Trans,1.,J,D3,0.,Tmp);
    blas::mm(NoTrans,NoTrans,1.,Tmp,J,0.,D4);

    GeMatrix<FullStorage<T,ColMajor> > m1(d_,d_);
    for (int j=1; j<=d_; ++j) {
        for (int k=1; k<=d_; ++k) {
            m1(j,k) = M1(d-2+j,d-2+k);
        }
    }

    GeMatrix<FullStorage<T,ColMajor> > m2(2*d_+d-3,d_);
    for (int j=1; j<=2*d_+d-3; ++j) {
        for (int k=1; k<=d_; ++k) {
            m2(j,k) = M1(d+d_-2+j,d-2+k);
        }
    }
    blas::mm(NoTrans,NoTrans,1.,D4,m1,0.,Tmp);

    blas::mm(NoTrans,NoTrans,1.,Tmp,inv(D4),0.,m1);
    Tmp = m2;
    blas::mm(NoTrans,NoTrans,1.,Tmp,inv(D4),0.,m2);

    for (int j=1; j<=d_; ++j) {
        for (int k=1; k<=d_; ++k) {
            MD(d-2+j,d-2+k) = m1(j,k);
        }
    }
    for (int j=1; j<=2*d_+d-3; ++j) {
        for (int k=1; k<=d_; ++k) {
            MD(d+d_-2+j,d-2+k) = m2(j,k);
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScaleDual_2(int d, int d_)
{
    int t = std::max(d+d_-2,4*d-6);
    GeMatrix<FullStorage<T,ColMajor> > TL(t,t);
    for (int j=1; j<=t; ++j) {
        TL(j,j) = 1.;
    }
    _twoScalePrimal(d,d_);
    _twoScaleDual_1(d,d_);
    MDD.engine().resize(3*d_+2*d-5,d+d_-2) || MDD.engine().fill();

    for (int j=1; j<=3*d_+2*d-5; ++j) {
        for (int k=1; k<=d+d_-2; ++k) {
            MDD(j,k) = MD(j,k);
        }
    }

    for (int k=d-2; k>=1; --k) {
        GeMatrix<FullStorage<T,ColMajor> > B(d-1+k,d-1+k), C;
        _orthonormalize(MP,MDD);
        for (int j=1; j<=d+d_-2; ++j) {
            for (int s=1; s<=d+d_-2; ++s) {
                TL(j,s) = TT(j,s);
            }
        }
        for (int j=1; j<=d-1+k; ++j) {
            for (int s=1; s<=d-1+k; ++s) {
                for (int l=1; l<=2*k+2*d-2; ++l) {
                    B(j,s) += MP(l,j)*TL(l,k+s);
                }
                B(j,s) *= .5;
            }
        }

        C = inv(B);

        MDD(k,k) = 0;
        DenseVector<Array<T> > x(d+k-1), y(d+k-1);
        for (int j=1; j<=d+k-1; ++j) {
            x(j) = MDD(k,k)*.5*MP(k,j);
        }
        y(k) = 1; x *= -1; x += y; y = C*x;
        for (int j=1; j<=d+k-1; ++j) {
            MDD(k+j,k) = y(j);
        }
    }
    _orthonormalize(MP,MDD);
    GeMatrix<FullStorage<T,ColMajor> > InvT, D(2*d_+d-3,2*d_+d-3), Tmp, TmpTT;
    InvT = inv(TT);

    GeMatrix<FullStorage<T,ColMajor> > K(3*d_+2*d-5,3*d_+2*d-5);
    K.diag(0) = 1;
    K(TT.rows(),TT.cols()) = TT;
    blas::mm(NoTrans,NoTrans,1.,K,MDD,0.,Tmp);
    blas::mm(NoTrans,NoTrans,1.,Tmp,InvT,0.,MDDD);
}


template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_calcM0_()
{
    _twoScaleDual_2(d,d_);

    GeMatrix<FullStorage<T,ColMajor> > Mj0_(pow2i<T>(min_j0+1)+d-1,
                                            pow2i<T>(min_j0)+d-1),
                                       Mj0_Right, Bugfix;
    T factor = Const<T>::R_SQRT2;
    Bugfix = MDDD;
    blas::scal(factor,Bugfix);
    Mj0_(MDDD.rows(), MDDD.cols()) = Bugfix;
    int row = d+d_-1;
    for (int c=MDDD.lastCol()+1; c<=Mj0_.lastCol()-MDDD.numCols(); ++c,row+=2) {
        Mj0_(_(row,row+phi_R.a_.length()-1),c) = factor*phi_R.a_;
    }
    arrow(Bugfix,Mj0_Right);
    Mj0_(_(Mj0_.lastRow()-Mj0_Right.numRows()+1,Mj0_.lastRow()),
         _(Mj0_.lastCol()-Mj0_Right.numCols()+1,Mj0_.lastCol())) = Mj0_Right;

    M0_ = RefinementMatrix<T,Interval,Primbs>(Mj0_Right.numCols(),Mj0_Right.numCols(), 
                                              Mj0_, min_j0, min_j0);
    setLevel(_j);
}

} // namespace lawa

