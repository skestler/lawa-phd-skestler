namespace lawa {



template <typename T>
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::AdaptivePDEOperatorOptimized1D
                                             (const ReallineCDFBasis1D &_basis, T _reaction,
                                              T _convection, T _diffusion, T /*thresh*/,
                                              int /*NumOfCols*/, int /*NumOfRows*/)
: basis(_basis), reaction(_reaction), convection(_convection),
  diffusion(_diffusion), cA(0.), CA(0.), kappa(0.),
  compression(basis), pde_op1d(basis,reaction,convection,diffusion), prec1d(),
  pde_data1d(pde_op1d, prec1d, compression),
  P_data()
{

}

template <typename T>
T
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::operator()(const Index1D &row_index,
                                                        const Index1D &col_index)
{
    T val = pde_data1d(row_index,col_index);
    return this->prec(row_index) * val * this->prec(col_index);
}

template <typename T>
T
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = pde_data1d(index,index);
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::mv(const IndexSet<Index1D> &LambdaRow,
                                                const Coefficients<Lexicographical,T,Index1D> &v)
{

    compression.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index1D> w_;
    for (const_coeff1d_it mu = v.begin(); mu != v.end(); ++mu) {
        Index1D col_index = (*mu).first;
        T prec_colindex = this->prec(col_index);
        IndexSet<Index1D> LambdaRowSparse = compression.SparsityPattern(col_index, LambdaRow);
        for (const_set1d_it lambda=LambdaRowSparse.begin(); lambda!=LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += pde_data1d((*lambda),col_index) * prec_colindex * (*mu).second;
        }
    }
    for (coeff1d_it it=w_.begin(); it!=w_.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }
    return w_;

    return lawa::mv_sparse(LambdaRow, (*this), v);
}

template <typename T>
void
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow,
                                                                 const IndexSet<Index1D>& LambdaCol,
                                                                 SparseMatrixT &A_flens, int J)
{
    //lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);


        std::map<Index1D,int,lt<Lexicographical,Index1D> > row_indices;
        int row_count = 1, col_count = 1;
        for (const_set1d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
        }
        this->compression.setParameters(LambdaRow);
        for (const_set1d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
            IndexSet<Index1D> LambdaRowSparse = this->compression.SparsityPattern(*col, LambdaRow, J);
            for (const_set1d_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
                T val = this->operator()(*row,*col);
                if (fabs(val)>0) {
                    A_flens(row_indices[*row],col_count) = val;
                }
            }
        }
        A_flens.finalize();
}

template <typename T>
void
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::toFlensSparseMatrix
                                                      (const IndexSet<Index1D>& LambdaRow,
                                                       const IndexSet<Index1D>& LambdaCol,
                                                       SparseMatrixT &A_flens, T eps)
{
    int J=0;        //compression
    T CA_temp = 5;
    J = std::min((int)std::ceil(-1./(basis.d-1.5)*log(eps/CA_temp)/log(2.)), 25);
    //J = std::min((int)std::ceil(-1./(basis.d-1.5)*log(eps/CA)/log(2.)), 25);
    std::cerr << "AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>: Estimated compression level for "
              << "tol = " << eps << " : " << J << std::endl;
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,J);
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                   int k, int J)
{
    std::cerr << "AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>: Constants cA and CA need to be "
              << "estimated first. Exit." << std::endl;

    int d=basis.d;
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index1D> temp;
    temp = v;

    int s = 0, count = 0;
    for (const_abs_coeff1d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        Index1D colindex = (*it).second;
        T prec_colindex = this->prec(colindex);
        IndexSet<Index1D> Lambda_v;
        int maxlevel;
        J==-1000 ? maxlevel=colindex.j+(k-s)+1 : maxlevel=J;
        Lambda_v=lambdaTilde1d_PDE(colindex, basis,(k-s), basis.j0, std::min(maxlevel,36), false);
        for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
            ret[*mu] += this->pde_data1d(*mu, (*it).second) * prec_colindex * (*it).first;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }
    for (coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    return ret;
}

template <typename T>
void
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                            T eps, Coefficients<Lexicographical,T,Index1D> &ret)
{
    std::cerr << "AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>: Constants cA and CA need to be "
              << "estimated first. Exit." << std::endl;


    //Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size()==0) return;// ret;

    // Implementation of APPLY as stated in Urban:2009, p. 216
/*
    Coefficients<AbsoluteValue,T,Index1D> v_abs;
    v_abs = v;
    int k = this->findK(v_abs, eps);
    ret = this->apply(v, k);
    return;// ret;
*/


    // Implementation of APPLY as stated in Stevenson:2009, p.560.

    Coefficients<Bucket,T,Index1D> v_bucket;
    T tol = 0.5*eps/CA;
    v_bucket.bucketsort(v,tol);
    long double squared_v_norm = (long double)std::pow(v.norm(2.),2.);
    long double squared_v_bucket_norm = 0.;
    T delta=0.;
    int l=0;
    int support_size_all_buckets=0;

    for (int i=0; i<(int)v_bucket.buckets.size(); ++i) {
        squared_v_bucket_norm += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
        T squared_delta = fabs(squared_v_norm - squared_v_bucket_norm);
        support_size_all_buckets += v_bucket.buckets[i].size();
        delta = std::sqrt(squared_delta);
        l = i+1;
        if (squared_delta<tol*tol) {
            break;
        }
    }
    // When tol*tol is close to machine precision, the above approach does not work any longer.
    // Also summing first over the smallest terms does not help (an implementation can be found at
    // the end of the file) -> i<=l required in the loop below!!!.
    // Coarsening does not suffer form this effect so that we have delta < eps/2.
    if (delta>tol) delta = eps/2.;

    for (int i=0; i<l; ++i) {
        Coefficients<Lexicographical,T,Index1D> w_p;
        v_bucket.addBucketToCoefficients(w_p,i);
        if (w_p.size()==0) continue;
        T numerator = w_p.norm(2.) * support_size_all_buckets;
        //T denominator = w_p.size() * 0.5*tol / CA;
        T denominator = w_p.size() * (eps-delta) / CA;
        //std::cout << "Bucket " << i << ": numerator=" << numerator << ", denominator=" << denominator << std::endl;
        if (denominator < 0) {
            //std::cout << "Bucket " << i << ": eps=" << eps << ", delta=" << delta << std::endl;
        }
        int jp = (int)std::max(std::log(numerator/denominator) / std::log(2.) / (basis.d-1.5), 0.);
        //std::cout << "Bucket " << i << ": #wp= " << w_p.size() << ", jp=" << jp << std::endl;
        jp = std::min(jp,36);
        for (const_coeff1d_it it=w_p.begin(); it!=w_p.end(); ++it) {
            Index1D colindex = (*it).first;
            T prec_colindex = this->prec(colindex);
            IndexSet<Index1D> Lambda_v;

            Lambda_v=lambdaTilde1d_PDE(colindex, basis,jp, basis.j0, std::min(colindex.j+jp,36), false);
            for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += this->pde_data1d(*mu, colindex) * prec_colindex * (*it).second;
            }


        }
    }
    for (coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    return;

}

template <typename T>
void
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                            T eps, const IndexSet<Index1D> &Lambda,
                                                            Coefficients<Lexicographical,T,Index1D> &ret)
{
    //todo: optimize!!
    this->apply(v,eps,ret);
    ret = P(ret,Lambda);
}

template <typename T>
int
AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>::findK(const Coefficients<AbsoluteValue,T,Index1D> &v,
                                                   T eps)
{
    int d=basis.d;
    if (v.size() == 0) return 1;
    T s=d-1.5;    //s = gamma-1, gamma the smoothness index of the wavelet basis

    T tau = 1.0 / (s + 0.5);
    // here the constant in (7.27) (KU-Wavelet) is estimated with 10
    int k_eps = static_cast<int>(10*log(std::pow(eps, -1.0/s)*std::pow(v.wtauNorm(tau), 1.0/s)) / log(2.0));
    DenseVector<Array<T> > normsec = v.norm_sections();
    T ErrorEstimateFactor = 1.;
    //std::cout << "eps = " << eps << ", k_eps = " << k_eps << std::endl;

    for (int k=1; k<=k_eps; ++k) {
        //std::cout << "At k = " << setw(3) << k;

        T R_k = 0.0;
        for (int i=k; i<=normsec.lastIndex()-1; ++i) {
            R_k += normsec(i+1);
        }
        R_k *= this->CA;
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k += std::pow(2.,-k*s) * normsec(1);
        //std::cout << ", R_k = " << setw(10) << R_k;

        for (int l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()-1) {
                //R_k += std::pow(l,-1.01)*std::pow(2.,-l*s) * normsec(k-l+1);
                R_k += std::pow(2.,-l*s) * normsec(k-l+1);
            }
        }
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k *= ErrorEstimateFactor;
        //std::cout << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10) << eps << endl;

        if (R_k<=eps) {
            std::cout << "   findK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            int maxlevel=22;
            if (d==2)         {    maxlevel=28; }
            else if (d==3)    {   maxlevel=21; }    //for non-singular examples, also lower values are possible
                                                    //high level for ex. 3, j0=-inf required.
            return std::min(k,maxlevel);
        }
    }
    return std::min(k_eps,25);    //higher level differences result in translation indices that cannot be stored in int.
}





template<typename T, DomainType Domain>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::AdaptivePDEOperatorOptimized1D
                                                           (const SparseMultiBasis1D &_basis,
                                                            T _reaction, T _convection, T _diffusion,
                                                            T /*thresh*/, int /*NumOfCols*/,
                                                            int /*NumOfRows*/)
: basis(_basis), reaction(_reaction), convection(_convection), diffusion(_diffusion),
  cA(0.), CA(0.), kappa(0.),
  compression(basis), pde_op1d(basis,reaction,convection,diffusion), prec1d(),
  pde_data1d(pde_op1d, prec1d, compression),
  P_data()
{
    std::cout << "AdaptivePDEOperatorOptimized1D<T, Primal, Domain, SparseMulti>: j0="
              << basis.j0 << std::endl;
    if (basis.d==4) {

        if      (basis.j0==0)  {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-1) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-2) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-3) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-4) {    cA = 0.25;  CA = 2.9;    }
        else assert(0);

    }
    kappa = CA/cA;
}

template<typename T, DomainType Domain>
T
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::operator()(const Index1D &row_index,
                                                                         const Index1D &col_index)
{
    if ( abs(row_index.j-col_index.j)>1) {
        return 0.;
    }
    else {
        T val = pde_data1d(row_index,col_index);
        return this->prec(row_index) * val * this->prec(col_index);
    }
}

template<typename T, DomainType Domain>
T
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = pde_data1d(index,index);
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::mv
                                                  (const IndexSet<Index1D> &LambdaRow,
                                                   const Coefficients<Lexicographical,T,Index1D> &v)
{
    Coefficients<Lexicographical,T,Index1D> ret;
    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        Index1D col_index = (*it).first;
        IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
        for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
            Index1D row_index = *set_it;
            if (LambdaRow.count(row_index)>0) {
                ret[*set_it] += this->operator()(row_index,col_index) * (*it).second;
            }
        }
    }

    return ret;
}

template<typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::toFlensSparseMatrix
                                                           (const IndexSet<Index1D>& LambdaRow,
                                                            const IndexSet<Index1D>& LambdaCol,
                                                            SparseMatrixT &A_flens, int /*J*/)
{
    std::map<Index1D,int,lt<Lexicographical,Index1D> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set1d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    this->compression.setParameters(LambdaRow);
    for (const_set1d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        IndexSet<Index1D> LambdaRowSparse = this->compression.SparsityPattern(*col, LambdaRow);
        for (const_set1d_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
            T val = this->operator()(*row,*col);
            if (fabs(val)>0) {
                A_flens(row_indices[*row],col_count) = val;
            }
        }
    }
    A_flens.finalize();
}

template <typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T,Primal,Domain,SparseMulti>::toFlensSparseMatrix
                                                        (const IndexSet<Index1D>& LambdaRow,
                                                         const IndexSet<Index1D>& LambdaCol,
                                                         SparseMatrixT &A_flens, T /*eps*/)
{
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,1);
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  int/* k */, int /* J */, cxxblas::Transpose trans)
{
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    if (trans==cxxblas::NoTrans) {
        for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
            Index1D col_index = (*it).first;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
            T prec_col_index = this->prec(col_index);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                Index1D row_index = *set_it;
                T val = pde_data1d(row_index,col_index) * prec_col_index * (*it).second;
                if (fabs(val)>0) {
                    ret[row_index] += val;
                }
            }
        }
    }
    else if (trans==cxxblas::Trans) {
        for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
            Index1D row_index = (*it).first;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(row_index, basis, 1, basis.j0);
            T prec_row_index = this->prec(row_index);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                Index1D col_index = *set_it;
                T val = pde_data1d(row_index,col_index) * prec_row_index * (*it).second;
                if (fabs(val)>0) {
                    ret[col_index] += val;
                }
            }
        }
    }
    for (const_coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= this->prec((*it).first);
    }
    return ret;
}

template<typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  T /*eps*/, Coefficients<Lexicographical,T,Index1D> &ret,
                                                  cxxblas::Transpose trans)
{
    //Coefficients<Lexicographical,T,Index1D> ret;
    ret = this->apply(v,0,0,trans);
    return;// ret;
}

template <typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T,Primal,Domain,SparseMulti>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                                   T eps, const IndexSet<Index1D> &Lambda,
                                                                   Coefficients<Lexicographical,T,Index1D> &ret,
                                                                   cxxblas::Transpose trans)
{
    //todo: optimize!!
    this->apply(v,eps,ret,trans);
    ret = P(ret,Lambda);
}




}
