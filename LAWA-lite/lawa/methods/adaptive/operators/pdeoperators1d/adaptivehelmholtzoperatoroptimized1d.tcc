namespace lawa {

template <typename T>
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::AdaptiveHelmholtzOperatorOptimized1D
                                             (const ReallineCDFBasis1D &_basis, bool _w_XBSpline,
                                              T _c, T /*thresh*/, int /*NumOfCols*/, int /*NumOfRows*/)
: basis(_basis), w_XBSpline(_w_XBSpline), c(_c),
  cA(0.), CA(0.), kappa(0.),
  compression(basis), helmholtz_op1d(basis,c), prec1d(),
  helmholtz_data1d(helmholtz_op1d, prec1d, compression),
  P_data()
{
    std::cerr << "AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>: j0=";
    if (!w_XBSpline) {    //only wavelet discretization
        std::cerr << "-Infinity" << std::endl;
        if (basis.d==2 && basis.d_==2 && c == 1.) {
            cA = 0.19; CA = 2.1;
        }
        else if (basis.d==3 && basis.d_==3 && c == 1.) {
            cA = 0.038; CA = 2.7;
        }
        else if (basis.d==3 && basis.d_==5 && c == 1.) {
            cA = 0.16; CA = 2.8;
        }
        else assert(0);
    }
    else {                //wavelet with b-splines on coarse level discretization
        std::cerr << basis.j0 << std::endl;
        if (basis.d==2 && basis.d_==2 && c == 1.) {
            if      (basis.j0==2)  {    cA = 0.375; CA = 2.1;    }
            else if (basis.j0==1)  {    cA = 0.375; CA = 2.1;    }
            else if (basis.j0==0)  {    cA = 0.375; CA = 2.1;    }
            else if (basis.j0==-1) {    cA = 0.58;  CA = 2.1;    }
            else if (basis.j0==-2) {    cA = 0.58;  CA = 2.1;    }
            else if (basis.j0==-3) {    cA = 0.55;  CA = 2.1;    }
            else if (basis.j0==-4) {    cA = 0.46;  CA = 2.1;    }
            else if (basis.j0==-5) {    cA = 0.4;   CA = 2.1;    }
            else if (basis.j0==-6) {    cA = 0.36;  CA = 2.1;    }
            else if (basis.j0==-7) {    cA = 0.33;  CA = 2.1;    }
            else if (basis.j0==-8) {    cA = 0.3;   CA = 2.1;    }
            else if (basis.j0==-8) {    cA = 0.29;  CA = 2.1;    }
            else if (basis.j0==-10) {   cA = 0.27;  CA = 2.1;    }
            else assert(0);
        }
        else if (basis.d==3 && basis.d_==3 && c == 1.) {
            if (basis.j0==0)  {    cA = 0.43;  CA = 1.94;    }
            else if (basis.j0==-1) {    cA = 0.39;  CA = 2.03;    }
            else if (basis.j0==-2) {    cA = 0.29;  CA = 2.24;    }
            else if (basis.j0==-3) {    cA = 0.21;  CA = 2.44;    }
            else if (basis.j0==-4) {    cA = 0.15;  CA = 2.55;    }
            else if (basis.j0==-5) {    cA = 0.12;  CA = 2.59;    }
            else if (basis.j0==-6) {    cA = 0.11;   CA = 2.61;    }
            else if (basis.j0==-7) {    cA = 0.09;  CA = 2.62;    }
            else if (basis.j0==-8) {    cA = 0.078; CA = 2.63;    }
            else if (basis.j0==-9) {    cA = 0.07;  CA = 2.64;    }
            else if (basis.j0==-10) {    cA = 0.063; CA = 2.64;    }
            else assert(0);
        }
        else if (basis.d==3 && basis.d_==5 && c == 1.) {
            if (basis.j0==0)  {    cA = 0.457;  CA = 1.96;    }
            else if (basis.j0==-1) {    cA = 0.416;  CA = 2.07;    }
            else if (basis.j0==-2) {    cA = 0.326;  CA = 2.32;    }
            else if (basis.j0==-3) {    cA = 0.24;   CA = 2.53;    }
            else if (basis.j0==-4) {    cA = 0.197;  CA = 2.66;    }
            else if (basis.j0==-5) {    cA = 0.18;   CA = 2.7;    }
            else if (basis.j0==-6) {    cA = 0.17;   CA = 2.71;    }
            else assert(0);
        }
        else assert(0);
    }
    kappa = CA/cA;
}

template <typename T>
T
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::operator()(const Index1D &row_index,
                                                        const Index1D &col_index)
{
    T val = helmholtz_data1d(row_index,col_index);
    return this->prec(row_index) * val * this->prec(col_index);
}

template <typename T>
T
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = helmholtz_data1d(index,index);
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::mv(const IndexSet<Index1D> &LambdaRow,
                                                const Coefficients<Lexicographical,T,Index1D> &v)
{

    compression.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index1D> w_;
    for (const_coeff1d_it mu = v.begin(); mu != v.end(); ++mu) {
        Index1D col_index = (*mu).first;
        T prec_colindex = this->prec(col_index);
        IndexSet<Index1D> LambdaRowSparse = compression.SparsityPattern(col_index, LambdaRow);
        for (const_set1d_it lambda=LambdaRowSparse.begin(); lambda!=LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += helmholtz_data1d((*lambda),col_index) * prec_colindex * (*mu).second;
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
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow,
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
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::toFlensSparseMatrix
                                                      (const IndexSet<Index1D>& LambdaRow,
                                                       const IndexSet<Index1D>& LambdaCol,
                                                       SparseMatrixT &A_flens, T eps)
{
    int J=0;        //compression
    J = std::min((int)std::ceil(-1./(basis.d-1.5)*log(eps/CA)/log(2.)), 25);
    std::cerr << "AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>: Estimated compression level for "
              << "tol = " << eps << " : " << J << std::endl;
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,J);
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                   int k, int J)
{
    int d=basis.d;
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index1D> temp;
    temp = v;

    if (w_XBSpline) {
        int s = 0, count = 0;
        for (const_abs_coeff1d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            Index1D colindex = (*it).second;
            T prec_colindex = this->prec(colindex);
            IndexSet<Index1D> Lambda_v;
            int maxlevel;
            J==-1000 ? maxlevel=colindex.j+(k-s)+1 : maxlevel=J;
            Lambda_v=lambdaTilde1d_PDE(colindex, basis,(k-s), basis.j0, std::min(maxlevel,36), false);
            for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += this->helmholtz_data1d(*mu, (*it).second) * prec_colindex * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
        for (coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
            (*it).second *= this->prec((*it).first);
        }
    }
    else {
        int s = 0, count = 0;
        for (const_abs_coeff1d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index1D> Lambda_v;
            Index1D colindex = (*it).second;
            int maxlevel=colindex.j+(k-s)+1;
            Lambda_v=lambdaTilde1d_PDE_WO_XBSpline(colindex, basis, (k-s), (*it).second.j-(k-s),
                                                   std::max(36,maxlevel));
            for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += this->operator()(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
    }

    return ret;
}

template <typename T>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                            T eps, Coefficients<Lexicographical,T,Index1D> &ret)
{
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
            if (w_XBSpline) {
                Lambda_v=lambdaTilde1d_PDE(colindex, basis,jp, basis.j0, std::min(colindex.j+jp,36), false);
                for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                    ret[*mu] += this->helmholtz_data1d(*mu, colindex) * prec_colindex * (*it).second;
                }
            }
            else {
                Lambda_v=lambdaTilde1d_PDE_WO_XBSpline(colindex, basis,jp, colindex.j-jp, std::min(colindex.j+jp,36));
                for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                    ret[*mu] += this->helmholtz_data1d(*mu, colindex) * prec_colindex * (*it).second;
                }
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
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                            T eps, const IndexSet<Index1D> &Lambda,
                                                            Coefficients<Lexicographical,T,Index1D> &ret)
{
    //todo: optimize!!
    this->apply(v,eps,ret);
    ret = P(ret,Lambda);
}

template <typename T>
int
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>::findK(const Coefficients<AbsoluteValue,T,Index1D> &v,
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
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::AdaptiveHelmholtzOperatorOptimized1D
                                                           (const MWBasis1D &_basis,
                                                            T _c, T /*thresh*/, int /*NumOfCols*/,
                                                            int /*NumOfRows*/)
: basis(_basis), c(_c),
  cA(0.), CA(0.), kappa(0.),
  compression(basis), laplace_op1d(basis), prec1d(),
  laplace_data1d(laplace_op1d, prec1d, compression),
  P_data()
{
    std::cout << "AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>: j0="
              << basis.j0 << std::endl;
    if (Domain == R) {
        if (basis.d==2 && c==1.) {
            if (basis.j0==0)       {    cA = 0.14;  CA = 2.9;    }
            else if (basis.j0==-1) {    cA = 0.14;  CA = 2.9;    }
            else if (basis.j0==-2) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-3) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-4) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-5) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-6) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-7) {    cA = 0.19;  CA = 2.9;    }
            else if (basis.j0==-8) {    cA = 0.19;  CA = 2.9;    }
            else assert(0);
        }
        else if (basis.d==3 && c==1.) {
            if (basis.j0==0)       {    cA = 0.02;  CA = 2.15;    }
            else if (basis.j0==-1) {    cA = 0.077;  CA = 2.15;    }
            else if (basis.j0==-2) {    cA = 0.22;  CA = 2.15;    }
            else if (basis.j0==-3) {    cA = 0.36;  CA = 2.15;    }
            else if (basis.j0==-4) {    cA = 0.36;  CA = 2.15;    }
            else assert(0);
        }
        else {
            //assert(0);
        }
    }
    else if (Domain == Interval) {
        if (basis.d==2 && c==0.) {
            if (basis.j0==0)       {    cA = 0.18;    CA = 2.7;    }
            else assert(0);
        }
        else if (basis.d==3 && c==0.) {
            if (basis.j0==0)       {    cA = 0.306;   CA = 2.08;    }
            else assert(0);
        }
        else if (basis.d==4 && c==0.) {
            if (basis.j0==0)       {    cA = 0.329;   CA = 1.96;    }
            else assert(0);
        }
        else assert(0);
    }
    else {
        //assert(0);
    }

    kappa = CA/cA;
}

template<typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::operator()(const Index1D &row_index,
                                                                    const Index1D &col_index)
{
    T val = laplace_data1d(row_index,col_index);
    if (row_index.k==col_index.k && row_index.j==col_index.j && row_index.xtype==col_index.xtype) {
//    if (row_index.val==col_index.val) {
       val += this->c;
    }
    return this->prec(row_index) * val * this->prec(col_index);
}

template<typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = laplace_data1d(index,index);
        val += this->c;
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::mv(const IndexSet<Index1D> &LambdaRow,
                                                           const Coefficients<Lexicographical,T,Index1D> &v)
{
    return mv_sparse(LambdaRow, (*this), v);
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::toFlensSparseMatrix
                                                        (const IndexSet<Index1D>& LambdaRow,
                                                         const IndexSet<Index1D>& LambdaCol,
                                                         SparseMatrixT &A_flens, int J)
{
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

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::toFlensSparseMatrix
                                                        (const IndexSet<Index1D>& LambdaRow,
                                                         const IndexSet<Index1D>& LambdaCol,
                                                         SparseMatrixT &A_flens, T eps)
{
    int J=0;        //compression
    J = (int)std::ceil(-1./(basis.d-1.5)*log(eps/CA)/log(2.));
    std::cerr << "AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>: Estimated compression level for "
              << "tol = " << eps << " : " << J << std::endl;
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,J);
}

template <typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::apply
                                                  (const Coefficients<Lexicographical,T,Index1D> &v,
                                                   int k, int J)
{
    int d=basis.d;
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index1D> temp;
    temp = v;

    int s = 0, count = 0;
    for (const_abs_coeff1d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        ret[(*it).second] += c * (*it).first * this->prec((*it).second);
        Index1D colindex = (*it).second;
        T prec_colindex = this->prec(colindex);
        IndexSet<Index1D> Lambda_v;
        int maxlevel;
        J==-1000 ? maxlevel=colindex.j+(k-s)+1 : maxlevel=J;
        //JMAX is defined in index.h
        Lambda_v=lambdaTilde1d_PDE(colindex, basis,(k-s), basis.j0, std::min(maxlevel,JMAX), false);
        for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
            T val = this->laplace_data1d(*mu, (*it).second);
            val *= prec_colindex * (*it).first;

            ret[*mu] += val;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }
    for (coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    return ret;
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::apply
                                                  (const Coefficients<Lexicographical,T,Index1D> &v,
                                                   T eps, Coefficients<Lexicographical,T,Index1D> & ret)
{
/*
    Coefficients<AbsoluteValue,T,Index1D> v_abs;
    v_abs = v;
    int k = this->findK(v_abs, eps);
    //Coefficients<Lexicographical,T,Index1D> ret;
    ret = this->apply(v, k);
    return;// ret;
*/
    T comprConst = CA;
    T convFactor = std::min(basis.d-1.5,1.5);
    /* multiwavelets are only C^1 and not C^{d-2} and so their smoothness index is
       gamma = std::min(d,3)-1/2. Moreoever, the half order of the operator order here is t=1.
       We have \| A-A_j \| \leq comprConst * 2^{-convFactor*j}
       Replacing j -> j/convFactor yields the setting from Stevenson:2009
     */

    Coefficients<Bucket,T,Index1D> v_bucket;
    T tol = 0.5*eps/comprConst;
    v_bucket.bucketsort(v,tol);
    //std::cerr << "APPLY: NumOfBuckets=" << v_bucket.buckets.size() << std::endl;
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

    if (delta>eps) delta = eps/2.;  // due to rounding error (see also above for CDF wavelets)

    for (int i=0; i<l; ++i) {
        Coefficients<Lexicographical,T,Index1D> w_p;
        v_bucket.addBucketToCoefficients(w_p,i);
        if (w_p.size()==0) continue;
        T numerator = w_p.norm(2.) * support_size_all_buckets;
        T denominator = w_p.size() * (eps-delta) / comprConst;
        //T denominator = w_p.size() * 0.5*tol / CA;
        //std::cout << "Bucket " << i << ": size=" << w_p.size() << ", (tol-delta)" << fabs(tol-delta) << std::endl;
        int jp = (int)std::max((std::log(numerator/denominator) / std::log(2.)) / convFactor, 0.);
        //std::cout << "Bucket " << i << ": #wp= " << w_p.size() << ", jp=" << jp << std::endl;
        for (const_coeff1d_it it=w_p.begin(); it!=w_p.end(); ++it) {
            Index1D colindex = (*it).first;
            ret[colindex] += c * (*it).second * this->prec(colindex);
            T prec_colindex = this->prec(colindex);
            IndexSet<Index1D> Lambda_v;
            Lambda_v=lambdaTilde1d_PDE(colindex, basis,jp, basis.j0, std::min(colindex.j+jp,25), false);
            for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += this->laplace_data1d(*mu, colindex) * prec_colindex * (*it).second;
            }
        }
    }
    for (coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    return;

}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                                       T eps, const IndexSet<Index1D> &Lambda,
                                                                       Coefficients<Lexicographical,T,Index1D> &ret)
{
    //todo: optimize!!
    this->apply(v,eps,ret);
    ret = P(ret,Lambda);
}

template <typename T, DomainType Domain>
int
AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi>::findK
                                                    (const Coefficients<AbsoluteValue,T,Index1D> &v,
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
            else if (d==3)    {   maxlevel=19; }    //for non-singular examples, also lower values are possible
            return std::min(k,maxlevel);
        }
    }
    return std::min(k_eps,25);    //higher level differences result in translation indices that cannot be stored in int.
}



template<typename T, DomainType Domain>
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::AdaptiveHelmholtzOperatorOptimized1D
                                                           (const SparseMultiBasis1D &_basis, T _c)
: basis(_basis), c(_c),
  cA(0.), CA(0.), kappa(0.),
  compression(basis), helmholtz_op1d(basis,c), prec1d(),
  helmholtz_data1d(helmholtz_op1d, prec1d, compression),
  P_data()
{
    std::cout << "AdaptiveHelmholtzOperatorOptimized1D<T, Primal, Domain, SparseMulti>: j0="
              << basis.j0  << ", c = " << c<< std::endl;
    if (basis.d==4 && c==1.) {

        if      (basis.j0==0)  {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0== 1) {    cA = 0.25;  CA = 2.9;    }  //todo: remove, only for debug!!
        else if (basis.j0== 2) {    cA = 0.25;  CA = 2.9;    }  //todo: remove!!
        else if (basis.j0== 3) {    cA = 0.25;  CA = 2.9;    }  //todo: remove!!
        else if (basis.j0==-1) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-2) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-3) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-4) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-5) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-6) {    cA = 0.25;  CA = 2.9;    }
        else assert(0);

    }
    kappa = CA/cA;
}

template<typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::operator()(const Index1D &row_index,
                                                                      const Index1D &col_index)
{
    if ( abs(row_index.j-col_index.j)>1) {
        return 0.;
    }
    else {
        T val = helmholtz_data1d(row_index,col_index);
        return this->prec(row_index) * val * this->prec(col_index);
    }
}

template<typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = helmholtz_data1d(index,index);
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::mv
                                                  (const IndexSet<Index1D> &LambdaRow,
                                                   const Coefficients<Lexicographical,T,Index1D> &v)
{
    std::cerr << "mv called." << std::endl;
    Coefficients<Lexicographical,T,Index1D> ret;
    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        Index1D col_index = (*it).first;
        IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE<T>(col_index, basis, 1, basis.j0);
        for (const_set1d_it set_it=LambdaRow.begin(); set_it!=LambdaRow.end(); ++set_it) {
        //for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
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
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::toFlensSparseMatrix
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
        //for (const_set1d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row) {
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
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,Domain,SparseMulti>::toFlensSparseMatrix
                                                        (const IndexSet<Index1D>& LambdaRow,
                                                         const IndexSet<Index1D>& LambdaCol,
                                                         SparseMatrixT &A_flens, T /*eps*/)
{
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,1);
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  int/* k */, int /* J */)
{
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        Index1D col_index = (*it).first;
        IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE<T>(col_index, basis, 1, basis.j0);
        T prec_col_index = this->prec(col_index);
        for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
            Index1D row_index = *set_it;
            T val = helmholtz_data1d(row_index,col_index) * prec_col_index * (*it).second;
            if (fabs(val)>0) {
                ret[row_index] += val;
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
AdaptiveHelmholtzOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  T /*eps*/, Coefficients<Lexicographical,T,Index1D> &ret)
{
    //Coefficients<Lexicographical,T,Index1D> ret;
    ret = this->apply(v,0,0);
    return;// ret;
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized1D<T,Primal,Domain,SparseMulti>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                                         T eps, const IndexSet<Index1D> &Lambda,
                                                                         Coefficients<Lexicographical,T,Index1D> &ret)
{
    //todo: optimize!!
    this->apply(v,eps,ret);
    ret = P(ret,Lambda);
}



}   // namespace lawa


/*
 * std::cerr.precision(16);
    for (int i=(int)v_bucket.buckets.size()-1; i>=0; --i) {
        squared_v_bucket_norm += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
    }
    std::cerr << "Norms: " << squared_v_norm - tol*tol << " < " << squared_v_bucket_norm
              << " < " << squared_v_norm  << std::endl;
    long double diff = squared_v_bucket_norm - (squared_v_norm - tol*tol);
    std::cerr << "diff = " << diff << std::endl;
    long double tmp = 0.L;
    for (int i=(int)v_bucket.buckets.size()-1; i>=0; --i) {
        tmp += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
        if (tmp>diff) {
            l=i;
            break;
        }
    }
    std::cerr << "tmp = " << tmp << ", l = " << l << std::endl;
    tmp = 0.L;
    for (int i=l; i>=0; --i) {
        support_size_all_buckets += v_bucket.buckets[i].size();
        tmp += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
    }
    std::cerr << "tmp = " << tmp << ", squared_v_norm = " << squared_v_norm << std::endl;
    delta = std::sqrt(fabs(squared_v_norm - tmp));
    std::cerr << "Norms: delta = " << delta << ", tol = " << tol << std::endl;
    std::cerr.precision(8);
 */


/*
 * Experimental apply for no minimal level.
 *
        int s = 0, count = 0;
        T beta1 = 1./(d-0.5); //gamma-1 im nenner kuerzt sich raus!!
        T beta2 = beta1 - 0.1;

        for (const_abs_coeff1d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index1D> Lambda_v;
            int kj = k-s;
            int minlevel, maxlevel;
            int j=(*it).second.j;
            if (j>=0) {
                maxlevel = j+kj;
                minlevel = j-kj;
                if (minlevel<0) {
                    minlevel =  int(std::floor(beta2*j-beta1*kj));
                }
            }
            else {//j<0
                minlevel = int(std::floor(j-beta1*kj));
                maxlevel = int(std::ceil(j+beta1*kj));
                if (maxlevel>0) maxlevel = int(std::ceil((j+beta1*kj)/beta2));
            }

            Lambda_v=lambdaTilde1d_PDE_WO_XBSpline((*it).second, basis, kj, minlevel,
                                                    std::min(maxlevel,36));
            for (const_set1d_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += this->operator()(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
    }
*/
