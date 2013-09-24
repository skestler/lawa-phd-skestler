#include <algorithm>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Side,Domain,Cons> &basis) {
    IndexSet<Index1D> ret, tmp;
    typedef typename IndexSet<Index1D>::const_iterator const_it;
    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        tmp = C((*lambda),c,basis);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Side,Domain,Cons> &basis, const int Jmax) {
    IndexSet<Index1D> ret, tmp;
    typedef typename IndexSet<Index1D>::const_iterator const_it;
    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        tmp = C((*lambda),c,basis, Jmax);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis) {
    IndexSet<Index1D> ret;
    index_cone(lambda,c,basis,ret);
    return ret;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis, const int Jmax) {
    IndexSet<Index1D> ret;
    index_cone(lambda,c,basis,ret, Jmax);
    return ret;
}


// Security zone interval
template <typename T, Construction Cons>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,Interval,Cons> &basis,
           IndexSet<Index1D> &ret)
{
    using std::min;
    using std::max;

    int j=lambda.j, jP1, k=lambda.k;
    XType xtype = lambda.xtype;
    int kMin_mra = basis.mra.rangeI(j).firstIndex(), kMax_mra = basis.mra.rangeI(j).lastIndex();
    if (lambda.xtype==XBSpline) {
        jP1=j;
        ret.insert(Index1D(j,std::max(k-2,kMin_mra),xtype));
        ret.insert(Index1D(j,std::max(k-1,kMin_mra),xtype));
        ret.insert(Index1D(j,std::min(k+1,kMax_mra),xtype));
        ret.insert(Index1D(j,std::min(k+2,kMax_mra),xtype));
    } else {
        jP1=j+1;
    }
    Support<T> supp = basis.generator(xtype).support(j,k);

    T zLambda=0.5*(supp.l2+supp.l1);
    Support<T> contractedSupp(c*supp.l1 + (1-c)*zLambda, c*supp.l2 + (1-c)*zLambda);

    int kMin = basis.rangeJ(jP1).firstIndex(), kMax = basis.rangeJ(jP1).lastIndex();
    int kStart = std::min(std::max(iceil<int>(contractedSupp.l1 * pow2i<T>(jP1)), kMin), kMax);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kStart))>0));
    while ((kStart-1 >= kMin) &&
           (overlap(contractedSupp, basis.psi.support(jP1,std::max(kStart-1, kMin)))>0)) {
        --kStart;
    }
    int kEnd = std::max(std::min(ifloor(contractedSupp.l2 * pow2i<T>(jP1)), kMax), kMin);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kEnd))>0));
    while ((kEnd+1 <= kMax) &&
          (overlap(contractedSupp, basis.psi.support(jP1,std::min(kEnd+1, kMax)))>0)) {
        ++kEnd;
    }

    for (int k=kStart; k<=kEnd; ++k) {
        ret.insert(Index1D(jP1,k,XWavelet));
    }
}

template <typename T, Construction Cons>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,Interval,Cons> &basis,
           IndexSet<Index1D> &ret, const int Jmax)
{
    using std::min;
    using std::max;
    
    int j=lambda.j, jP1, k=lambda.k;
    XType xtype = lambda.xtype;
    int kMin_mra = basis.mra.rangeI(j).firstIndex(), kMax_mra = basis.mra.rangeI(j).lastIndex();
    if (lambda.xtype==XBSpline) {
        jP1=j;
        ret.insert(Index1D(j,std::max(k-2,kMin_mra),xtype));
        ret.insert(Index1D(j,std::max(k-1,kMin_mra),xtype));
        ret.insert(Index1D(j,std::min(k+1,kMax_mra),xtype));
        ret.insert(Index1D(j,std::min(k+2,kMax_mra),xtype));
    } else {
        jP1=j+1;
        if(jP1 >= Jmax){ 
            return;
        }
    }
    Support<T> supp = basis.generator(xtype).support(j,k);
    
    T zLambda=0.5*(supp.l2+supp.l1);
    Support<T> contractedSupp(c*supp.l1 + (1-c)*zLambda, c*supp.l2 + (1-c)*zLambda);
    
    int kMin = basis.rangeJ(jP1).firstIndex(), kMax = basis.rangeJ(jP1).lastIndex();
    int kStart = std::min(std::max(iceil<int>(contractedSupp.l1 * pow2i<T>(jP1)), kMin), kMax);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kStart))>0));
    while ((kStart-1 >= kMin) &&
           (overlap(contractedSupp, basis.psi.support(jP1,std::max(kStart-1, kMin)))>0)) {
        --kStart;
    }
    int kEnd = std::max(std::min(ifloor(contractedSupp.l2 * pow2i<T>(jP1)), kMax), kMin);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kEnd))>0));
    while ((kEnd+1 <= kMax) &&
           (overlap(contractedSupp, basis.psi.support(jP1,std::min(kEnd+1, kMax)))>0)) {
        ++kEnd;
    }
    
    for (int k=kStart; k<=kEnd; ++k) {
        ret.insert(Index1D(jP1,k,XWavelet));
    }
}

// Security zone periodic
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,Periodic,CDF> &basis,
           IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    //ret.insert(Index1D(j,k,xtype));
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,(k-1 >= basis.mra.rangeI(j).firstIndex()) ? k-1 : basis.mra.rangeI(j).lastIndex()
                                    + ((1 - (basis.mra.rangeI(j).firstIndex() - k+1))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+1 <= basis.mra.rangeI(j).lastIndex() ? k+1 : basis.mra.rangeI(j).firstIndex()
                                    - ((1 - (k+1 - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k-2 >= basis.mra.rangeI(j).firstIndex() ? k-2 : basis.mra.rangeI(j).lastIndex()
                                    + ((1 - (basis.mra.rangeI(j).firstIndex() - k+2))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+2 <= basis.mra.rangeI(j).lastIndex() ? k+2 : basis.mra.rangeI(j).firstIndex()
                                    - ((1 - (k+2 - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j)),xtype));
        Support<T> contractedSupp, supp = basis.mra.phi.phiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.psiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;
        /*    no wavelet indices on the same level?!
        long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
        */

        long int kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j+1,k1))>0)
            {
                int k = (int)k1;
                if(k < basis.rangeJ(j+1).firstIndex()){
                    k = basis.rangeJ(j+1).lastIndex() + ((1 - (basis.rangeJ(j+1).firstIndex() - k))%basis.cardJ(j+1));
                }
                if(k > basis.rangeJ(j+1).lastIndex()){
                    k = basis.rangeJ(j+1).firstIndex() - ((1 - (k - basis.rangeJ(j+1).lastIndex()))%basis.cardJ(j+1));
                }
                ret.insert(Index1D(j+1,k,XWavelet));
            }
        }

    }
}

template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,Periodic,CDF> &basis,
           IndexSet<Index1D> &ret, const int Jmax)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    //ret.insert(Index1D(j,k,xtype));
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,(k-1 >= basis.mra.rangeI(j).firstIndex()) ? k-1 : basis.mra.rangeI(j).lastIndex()
                           + ((1 - (basis.mra.rangeI(j).firstIndex() - k+1))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+1 <= basis.mra.rangeI(j).lastIndex() ? k+1 : basis.mra.rangeI(j).firstIndex()
                           - ((1 - (k+1 - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k-2 >= basis.mra.rangeI(j).firstIndex() ? k-2 : basis.mra.rangeI(j).lastIndex()
                           + ((1 - (basis.mra.rangeI(j).firstIndex() - k+2))%basis.mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+2 <= basis.mra.rangeI(j).lastIndex() ? k+2 : basis.mra.rangeI(j).firstIndex()
                           - ((1 - (k+2 - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j)),xtype));
        Support<T> contractedSupp, supp = basis.mra.phi.phiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;
        
        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);
        
        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    else {
        if(j+1 >= Jmax){ 
            return;
        }
        
        Support<T> contractedSupp, supp = basis.psi.psiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;
        /*    no wavelet indices on the same level?!
         long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
         long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);
         
         for (long int k1=kMin; k1<=kMax; ++k1) {
         if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
         int k = k1;
         if(k < basis.rangeJ(j).firstIndex()){
         k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
         }
         if(k > basis.rangeJ(j).lastIndex()){
         k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
         }
         ret.insert(Index1D(j,k,XWavelet));
         }
         }
         */
        
        long int kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);
        
        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j+1,k1))>0)
            {
                int k = (int)k1;
                if(k < basis.rangeJ(j+1).firstIndex()){
                    k = basis.rangeJ(j+1).lastIndex() + ((1 - (basis.rangeJ(j+1).firstIndex() - k))%basis.cardJ(j+1));
                }
                if(k > basis.rangeJ(j+1).lastIndex()){
                    k = basis.rangeJ(j+1).firstIndex() - ((1 - (k - basis.rangeJ(j+1).lastIndex()))%basis.cardJ(j+1));
                }
                ret.insert(Index1D(j+1,k,XWavelet));
            }
        }
        
    }
}

// Security zone realline
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis, IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,k-1,xtype));
        ret.insert(Index1D(j,k+1,xtype));
        ret.insert(Index1D(j,k-2,xtype));
        ret.insert(Index1D(j,k+2,xtype));

        Support<T> contractedSupp, supp = basis.mra.phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long int kMin, kMax;

        kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j+1,k1))>0) ret.insert(Index1D(j+1,k1,XWavelet));
        }
    }
}

// Security zone orthonormal multi-wavelet realline
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Orthogonal,R,Multi> &basis,
           IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;

    const BSpline<T,Orthogonal,R,Multi> &phi = basis.mra.phi;
    const Wavelet<T,Orthogonal,R,Multi> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;
    Support<T> max_support_refbspline = phi.max_support();
    Support<T> max_support_refwavelet = psi.max_support();

    if (xtype==XBSpline) {
        for (int i=1; i<=2*numSplines; ++i) {
            ret.insert(Index1D(j,k-i,xtype));
            ret.insert(Index1D(j,k+i,xtype));
        }

        Support<T> contractedSupp, supp = phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - psi.max_support().l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - psi.max_support().l1);
        for (int k_help=kMin; k_help<=kMax; ++k_help) {
            for (int k1=(k_help-1)*numWavelets+1; k1<=k_help*numWavelets; ++k1) {
                if (overlap(contractedSupp, basis.psi.support(j,k1))>0) {
                    ret.insert(Index1D(j,k1,XWavelet));
                }
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long int kMin, kMax;

        kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - psi.max_support().l2);
        kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - psi.max_support().l1);
        for (long int k_help=kMin; k_help<=kMax; ++k_help) {
            for (long int k1=(k_help-1)*numWavelets+1; k1<=k_help*numWavelets; ++k1) {
                if (overlap(contractedSupp, psi.support(j+1,k1))>0) {
                    ret.insert(Index1D(j+1,k1,XWavelet));
                }
            }
        }
    }
}

// Security zone orthonormal multi-wavelet interval
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Orthogonal,Interval,Multi> &basis,
           IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;

    const BSpline<T,Orthogonal,Interval,Multi> &phi = basis.mra.phi;
    const Wavelet<T,Orthogonal,Interval,Multi> &psi = basis.psi;

    if (xtype==XBSpline) {
        Support<T> supp = phi.support(j,k);
        int j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.getScalingNeighborsForScaling(j, k, basis,
                                            j_scaling, k_scaling_first, k_scaling_last);

        for (int k_scaling = k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            if (overlap(supp, phi.support(j_scaling,k_scaling))>0) {
                ret.insert(Index1D(j_scaling,k_scaling,XBSpline));
            }
        }

        int j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.getWaveletNeighborsForScaling(j, k, basis,
                                            j_wavelet, k_wavelet_first, k_wavelet_last);
        assert(j_wavelet==j);
        k_wavelet_first = std::max((long)basis.rangeJ(j_wavelet).firstIndex(), k_wavelet_first-basis._numInnerParts);
        k_wavelet_last  = std::min((long)basis.rangeJ(j_wavelet).lastIndex(),  k_wavelet_last+basis._numInnerParts);
        for (int k_wavelet = k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            if (overlap(supp, psi.support(j_wavelet,k_wavelet))>0) {
                ret.insert(Index1D(j_wavelet,k_wavelet,XWavelet));
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        int j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.getHigherWaveletNeighborsForWavelet(j, k, basis,
                                                  j_wavelet, k_wavelet_first, k_wavelet_last);
        assert(j_wavelet==j+1);
        k_wavelet_first = std::max((long)basis.rangeJ(j_wavelet).firstIndex(), k_wavelet_first-basis._numInnerParts);
        k_wavelet_last  = std::min((long)basis.rangeJ(j_wavelet).lastIndex(),  k_wavelet_last+basis._numInnerParts);
        for (int k_wavelet = k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            if (j_wavelet >= 28) {
                std::cerr << "Level too large!!" << std::endl;
            }
            if (overlap(supp, psi.support(j_wavelet,k_wavelet))>0) {
                ret.insert(Index1D(j_wavelet,k_wavelet,XWavelet));
            }
        }
    }
}

// Security zone special multi-wavelet realline
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,R,SparseMulti> &basis,
           IndexSet<Index1D> &ret)
{
    int  j=lambda.j;
    if (j>=JMAX) return;
    long  k=lambda.k;
    XType xtype=lambda.xtype;

    const BSpline<T,Primal,R,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,R,SparseMulti> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;
    Support<T> max_support_refwavelet = psi.max_support();

    //std::cout << "Index cone for " << lambda << std::endl;
    if (xtype==XBSpline) {
        for (int i=1; i<=2*numSplines; ++i) {
            ret.insert(Index1D(j,k-i,xtype));
            ret.insert(Index1D(j,k+i,xtype));
        }

        Support<T> contractedSupp, supp = phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        long kMin = floor( pow2i<long double>(j)*contractedSupp.l1 - max_support_refwavelet.l2) - 1;
        long kMax =  ceil( pow2i<long double>(j)*contractedSupp.l2 - max_support_refwavelet.l1) + 1;

        for (int k_row=kMin*numSplines; k_row<=kMax*numSplines; ++k_row) {
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << psi.support(j,k_row) << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long kMin, kMax;

        kMin = floor( pow2i<long double>(j+1)*contractedSupp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j+1)*contractedSupp.l2 - max_support_refwavelet.l1) / 2 + 1;

        kMin = kMin*numWavelets;
        kMax = kMax*numWavelets;

        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j+1,k_row);
            //std::cerr << "  -> wavelet (" << j+1 << ", " << k_row << "): " << supp_row << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, supp_row) > 0)  {
                ret.insert(Index1D(j+1,k_row,XWavelet));
            }
        }
        //std::cerr << "index cone: loop finished." << std::endl;
    }
}

// Security zone special multi-wavelet RPlus
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,RPlus,SparseMulti> &basis,
           IndexSet<Index1D> &ret)
{
    int  j=lambda.j;
    if (j>=JMAX) return;
    long  k=lambda.k;
    XType xtype=lambda.xtype;

    const BSpline<T,Primal,RPlus,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,RPlus,SparseMulti> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;
    Support<T> max_support_refwavelet = psi.max_support();

    //std::cout << "Index cone for " << lambda << std::endl;
    if (xtype==XBSpline) {
        long kMin = std::max(k-2*numSplines, basis.mra.rangeIL(j).firstIndex());
        for (int k_row=kMin; k_row<=k-1; ++k_row) {
            ret.insert(Index1D(j,k_row,xtype));
        }
        long kMax = k+2*numSplines;
        for (int k_row=k+1; k_row<=kMax; ++k_row) {
            ret.insert(Index1D(j,k_row,xtype));
        }

        Support<T> contractedSupp, supp = phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        kMin = floor( pow2i<long double>(j)*contractedSupp.l1 - max_support_refwavelet.l2) - 1;
        kMax =  ceil( pow2i<long double>(j)*contractedSupp.l2 - max_support_refwavelet.l1) + 1;
        kMin = std::max(kMin*numWavelets, basis.rangeJL(j).firstIndex());
        kMax = kMax*numWavelets;
        for (int k_row=kMin; k_row<=kMax; ++k_row) {
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << psi.support(j,k_row) << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long kMin, kMax;
        kMin = floor( pow2i<long double>(j+1)*contractedSupp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j+1)*contractedSupp.l2 - max_support_refwavelet.l1) / 2 + 1;
        kMin = std::max(kMin*numWavelets, basis.rangeJL(j).firstIndex());
        kMax = kMax*numWavelets;

        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j+1,k_row);
            //std::cerr << "  -> wavelet (" << j+1 << ", " << k_row << "): " << supp_row << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, supp_row) > 0)  {
                ret.insert(Index1D(j+1,k_row,XWavelet));
            }
        }
        //std::cerr << "index cone: loop finished." << std::endl;
    }
}

// Security zone special multi-wavelet interval
template <typename T>
void
index_cone(const Index1D &lambda, T c, const Basis<T,Primal,Interval,SparseMulti> &basis,
           IndexSet<Index1D> &ret)
{
    int  j=lambda.j;
    if (j>=JMAX) return;
    long  k=lambda.k;
    XType xtype=lambda.xtype;

    const BSpline<T,Primal,Interval,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,Interval,SparseMulti> &psi = basis.psi;
    int numSplines = (int)basis.mra._numSplines;
    int numWavelets = (int)basis._numSplines;
    Support<T> max_support_refwavelet = basis.max_support();

    //std::cout << "Index cone for " << lambda << std::endl;
    if (xtype==XBSpline) {
        long kMin = std::max(k-2*numSplines, basis.mra.long_rangeI(j).firstIndex());
        for (int k_row=kMin; k_row<=k-1; ++k_row) {
            ret.insert(Index1D(j,k_row,xtype));
        }
        long kMax = std::min(k+2*numSplines, basis.mra.long_rangeI(j).lastIndex());
        for (int k_row=k+1; k_row<=kMax; ++k_row) {
            ret.insert(Index1D(j,k_row,xtype));
        }

        Support<T> contractedSupp, supp = phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        kMin = floor( pow2i<long double>(j)*contractedSupp.l1 - max_support_refwavelet.l2) - 1;
        kMax =  ceil( pow2i<long double>(j)*contractedSupp.l2 - max_support_refwavelet.l1) + 1;
        kMin = std::max(kMin*numWavelets, basis.long_rangeJ(j+1).firstIndex());
        kMax = std::min(kMax*numWavelets, basis.long_rangeJ(j+1).lastIndex());
        for (int k_row=kMin; k_row<=kMax; ++k_row) {
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << psi.support(j,k_row) << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long kMin, kMax;
        kMin = floor( pow2i<long double>(j+1)*contractedSupp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j+1)*contractedSupp.l2 - max_support_refwavelet.l1) / 2 + 1;
        kMin = std::max(kMin*numWavelets, basis.long_rangeJ(j+1).firstIndex());
        kMax = std::min(kMax*numWavelets, basis.long_rangeJ(j+1).lastIndex());

        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j+1,k_row);
            //std::cerr << "  -> wavelet (" << j+1 << ", " << k_row << "): " << supp_row << " vs. " << contractedSupp << std::endl;
            if (overlap(contractedSupp, supp_row) > 0)  {
                ret.insert(Index1D(j+1,k_row,XWavelet));
            }
        }
        //std::cerr << "index cone: loop finished." << std::endl;
    }
}


template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,R,CDF> &basis, bool only_pos)
{
    IndexSet<Index1D> ret, tmp;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        if (only_pos && (*lambda).j<=0) continue;
        tmp = C_WO_XBSpline((*lambda),c,basis);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis) {
    const Wavelet<T,Primal,R,CDF> psi(basis,0);
    int j = lambda.j, k = lambda.k;
    IndexSet<Index1D> ret;
    Support<T> contractedSupp, supp = basis.psi.support(j,k);
    T center = 0.5*(supp.l1 + supp.l2);
    contractedSupp.l1 = c*supp.l1 + (1-c)*center;
    contractedSupp.l2 = c*supp.l2 + (1-c)*center;
    long int kMin, kMax;

    kMin = floor( pow2i<T>(j-1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j-1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j-1,k1))>0) ret.insert(Index1D(j-1,k1,XWavelet));
    }

    kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
    }

    kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j+1,k1))>0) ret.insert(Index1D(j+1,k1,XWavelet));
    }
    return ret;
}


//Security zone 2D
template <typename T, typename Basis2D>
IndexSet<Index2D>
C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis, bool extralevel)
{
    typedef typename IndexSet<Index2D>::const_iterator const_set2d_it;
    typedef typename IndexSet<Index1D>::const_iterator const_set1d_it;

    const_set2d_it Lambda_end = Lambda.end();
    IndexSet<Index2D>  ret;


    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_set2d_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*lambda).index1, (T)c, basis.first);
        C_index2 = C((*lambda).index2, (T)c, basis.second);

        if (!extralevel) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index2D new_index2d((*it_C_index1), (*lambda).index2);
                if (Lambda.find(new_index2d)==Lambda_end) {
                    ret.insert(new_index2d);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D new_index2d((*lambda).index1, (*it_C_index2));
                if (Lambda.find(new_index2d)==Lambda_end) {
                    ret.insert(new_index2d);
                }
            }
        }
        else {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    if (Lambda.count(Index2D((*it_C_index1), (*it_C_index2)))==0) {
                        ret.insert(Index2D((*it_C_index1), (*it_C_index2)));
                    }
                }
            }
        }
    }
    return ret;
}

//Security zone 2D with maximal level
template <typename T, typename Basis2D>
IndexSet<Index2D>
C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis, const int J1_max, const int J2_max )
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;
    
    IndexSet<Index2D>  ret;
    
    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*lambda).index1, c, basis.first, J1_max);
        C_index2 = C((*lambda).index2, c, basis.second, J2_max);
        
        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index2D((*lambda).index1, (*it_C_index2)))==0) {
                ret.insert(Index2D((*lambda).index1, (*it_C_index2)));
            }
        }
        /*
         for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
         for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
         if (Lambda.count(Index2D((*it_C_index1), (*it_C_index2)))==0) {
         ret.insert(Index2D((*it_C_index1), (*it_C_index2)));
         }
         }
         }
         */
    }
    return ret;
}

template <typename T, typename Basis2D>
IndexSet<Index2D>
C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis)
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    IndexSet<Index2D>  ret;

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*lambda).index1, c, basis.first);
        //C_index2 = C((*lambda).index2, c, basis.second);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
    }
    return ret;
}


//Security zone 3D
template <typename T, typename Basis3D>
IndexSet<Index3D>
C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis)
{
    typedef typename IndexSet<Index3D>::const_iterator const_it_3d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    IndexSet<Index3D>  ret;

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_3d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2, C_index3;
        C_index1 = C((*lambda).index1, c, basis.first);
        C_index2 = C((*lambda).index2, c, basis.second);
        C_index3 = C((*lambda).index3, c, basis.third);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index3D((*it_C_index1), (*lambda).index2, (*lambda).index3) )==0) {
                ret.insert(Index3D((*it_C_index1),   (*lambda).index2,  (*lambda).index3) );
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index3D((*lambda).index1, (*it_C_index2), (*lambda).index3) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*it_C_index2), (*lambda).index3) );
            }
        }
        for (const_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
            if (Lambda.count(Index3D((*lambda).index1, (*lambda).index2, (*it_C_index3) ) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*lambda).index2, (*it_C_index3) ) );
            }
        }
    }
    return ret;
}

} // namespace lawa

