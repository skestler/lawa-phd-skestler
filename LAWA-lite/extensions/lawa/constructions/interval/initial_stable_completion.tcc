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

#include <lawa/aux/densify.h>
#include <extensions/flens/lapack_flens.h>

namespace lawa {

template <typename T, Construction ConsPrimal, Construction ConsDual>
void
initial_stable_completion(const MRA<T,Primal,Interval,ConsPrimal> &mra,
                          const MRA<T,Dual,Interval,ConsDual> &mra_,
                          GeMatrix<FullStorage<T,cxxblas::ColMajor> > &M1,
                          GeMatrix<FullStorage<T,cxxblas::ColMajor> > &M1_)
{
    typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

    const RefinementMatrix<T,Interval,ConsPrimal> &M0 = mra.M0;
    const int d = mra.d;
    const int d_ = mra_.d_;
    const int j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
    const int l1 = mra.l1, l2 = mra.l2;
    M0.setLevel(j);

    
//    FullColMatrix *MM = new FullColMatrix[d/2];
    FullColMatrix MM[20];
    MRA<T,Primal,Interval,ConsPrimal> _mra(mra.d, j); // TODO BE ELIMINATED!!!!!!!!!!!!!!!!!!!!
    densify(cxxblas::NoTrans, _mra.M0, MM[0], 1,1);
    const int size = MM[0].numRows();
    
    FullColMatrix H(size,size);
    FullColMatrix HInv(size,size);
    FullColMatrix F(size, pow2i<T>(j));
    H.diag(0) = 1;
    HInv.diag(0) = 1;

    if (d&1) {
        for (int i=0; i<=(d-1)/2-1; ++i) {
            FullColMatrix Hi(size,size);
            Hi.diag(0) = 1;
            for (int p=1; p<=size; ++p) {
                for (int q=1; q<=size; ++q) {
                    if ((p==q+1) && ((p&1)==(i&1)) && (i<p) && (p<size-d+i+2)) {
                        Hi(p,q) = -MM[i](p,(p+i)/2) / MM[i](p-1,(p+i)/2);
                    }
                    if ((p==q-1) && ((p&1)==((i-d)&1)) && (d-i<=p) && (p<size-i)) {
                        Hi(p,q) = -MM[i](p,(p+d-i)/2) / MM[i](p+1,(p+d-i)/2);
                    }
                }
            }
            FullColMatrix HiInv(size,size);
            for (int p=1; p<=size; ++p) {
                for (int q=1; q<=size; ++q) {
                    if (p==q) {
                        if (((p&1)==(i&1)) && (p>1)) {
                            HiInv(p,q) = 1 / (1 - Hi(p-1,p)*Hi(p,p-1));
                        } else if (((p&1)!=(i&1)) && (p<pow2i<T>(j+1)+d-1-1)) {
                            HiInv(p,q) = 1 / (1 - Hi(p,p+1)*Hi(p+1,p));                            
                        } else {
                            HiInv(p,q) = 1;
                        }
                    }

                    if ((p==q+1) && ((p&1)==(i&1)) && (p>1)) {
                        HiInv(p,q) = -Hi(p,q) / (1 - Hi(p-1,p)*Hi(p,p-1));
                    }
                    if ((p==q-1) && ((p&1)==((i-d)&1)) && (p<(pow2i<T>(j+1)+d-1)-1)) {
                        HiInv(p,q) = -Hi(p,q) / (1 - Hi(p,p+1)*Hi(p+1,p));
                    }
                }
            }
            FullColMatrix HTmp = H, HInvTmp = HInv;
            blas::mm(NoTrans,NoTrans,1.,Hi,MM[i],0.,MM[i+1]);
            blas::mm(NoTrans,NoTrans,1.,Hi,HTmp,0.,H);
            blas::mm(NoTrans,NoTrans,1.,HInvTmp,HiInv,0.,HInv);
            Hi.engine().fill(0);
        }

        for (int p=F.firstRow(); p<=F.lastRow(); ++p) {
            for (int q=F.firstCol(); q<=F.lastCol(); ++q) {
                if ((p==2*q-1+(d-1)/2) && (q<=pow2i<T>(j-1))) {
                    F(p,q) = MM[(d-1)/2](p+1,q+(d-1)/2);
                }
                if ((p==2*q-1+(d+1)/2) && (q<=pow2i<T>(j-1))) {
                    F(p,q) = -MM[(d-1)/2](p-1,q+(d-1)/2);
                }
                if ((size-p+1==2*(pow2i<T>(j)-q+1)-1+(d-1)/2) 
                 && (q>pow2i<T>(j-1))) {
                    F(p,q) = MM[(d-1)/2](p-1,q+(d-1)/2);
                }
                if ((size-p+1==2*(pow2i<T>(j)-q+1)-1+(d+1)/2) 
                 && (q>pow2i<T>(j-1))) {
                    F(p,q) = -MM[(d-1)/2](p+1,q+(d-1)/2);
                }
            }
        }
    }

    if ((d&1)==0) {
        for (int i=0; i<=d/2-1; ++i) {
            FullColMatrix Hi(size,size);
            Hi.diag(0) = 1;
            for (int p=1; p<=size; ++p) {
                for (int q=1; q<=size; ++q) {
                    if ((p==q+1) && ((p&1)==(i&1)) && (i<p) && (p<=size-d+i+1)) {
                        Hi(p,q) = -MM[i](p,(p+i)/2) / MM[i](p-1,(p+i)/2);
                    }
                    if ((p==q-1) && ((p&1)==((i-d)&1)) && (d-1-i<p) && (p<=size-i)) {
                        Hi(p,q) = -MM[i](p,(p+d-i)/2) / MM[i](p+1,(p+d-i)/2);
                    }
                }
            }
            FullColMatrix HiInv(size,size);
            HiInv.diag(0) = 2.;
            cxxblas::geaxpy(cxxblas::ColMajor,NoTrans,size,size,-1.,
                            Hi.engine().data(),size,
                            HiInv.engine().data(),size);
            FullColMatrix HTmp = H, HInvTmp = HInv;
            blas::mm(NoTrans,NoTrans,1.,Hi,MM[0],0.,MM[i+1]);
            blas::mm(NoTrans,NoTrans,1.,Hi,HTmp,0.,H);
            blas::mm(NoTrans,NoTrans,1.,HInvTmp,HiInv,0.,HInv);
            Hi.engine().fill(0);
        }
        
        for (int p=F.firstRow(); p<=F.lastRow(); ++p) {
            for (int q=F.firstCol(); q<=F.lastCol(); ++q) {
                if (p==2*q-1+d/2) {
                    F(p,q) = 1.;
                }
            }
        }
    }
//    delete[] MM;
    
    //--- finally setting up M1 and M1_ ----------------------------------------
    //--- M1
    FullColMatrix Mj1Tmp, Mj1initial;
    blas::mm(NoTrans,NoTrans, 1., HInv, F, 0., Mj1Tmp);
    Mj1initial = Mj1Tmp(_(Mj1Tmp.firstRow()+mra_._bc(0), 
                          Mj1Tmp.lastRow()-mra_._bc(1)), _ );

    FullColMatrix Mj0, Mj0_;
    const RefinementMatrix<T,Interval,ConsDual> &M0_ = mra_.M0_;
    M0_.setLevel(j);
    densify(cxxblas::NoTrans, M0, Mj0, M0.firstRow(), M0.firstCol());
    densify(cxxblas::NoTrans, M0_, Mj0_, M0.firstRow(), M0.firstCol());
    FullColMatrix Tmp(M0.numRows(), M0_.numRows());
    Tmp.diag(0) = 1;
    blas::mm(NoTrans,Trans, -1., Mj0, Mj0_, 1., Tmp);
    blas::mm(NoTrans, NoTrans, 1., Tmp, Mj1initial, 0., M1);

    //--- M1_
    FullColMatrix Mj(pow2i<T>(j+1)+l2-l1-1-mra_._bc(0)-mra_._bc(1),
                     pow2i<T>(j+1)+l2-l1-1-mra_._bc(0)-mra_._bc(1));
    Mj( _ , _(1,M0.numCols())) = Mj0;
    Mj( _ , _(M0.numCols()+1, Mj.lastCol())) = M1;
    DenseVector<Array<int> > p(Mj.numRows());
    FullColMatrix MjInv, MjInvTmp = Mj, TransTmp2;
    
    MjInv = inv(Mj);

    FullColMatrix Mj1_Tmp = MjInv(_(pow2i<T>(j)+d-1+1-mra_._bc(0)-mra_._bc(1),
                                    pow2i<T>(j+1)+d-1-mra_._bc(0)-mra_._bc(1)), _ );
    copy(Trans, Mj1_Tmp, M1_);

    M1.engine().changeIndexBase(1+mra_._bc(0),1);
    M1_.engine().changeIndexBase(1+mra_._bc(0),1);

    // TODO: theoretical foundation for scaling this way!
    blas::scal(Const<T>::R_SQRT2,M1_);
    blas::scal(2*Const<T>::R_SQRT2,M1);
}

} // namespace lawa

