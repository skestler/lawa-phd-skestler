namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GHS_ADWAV(AdaptiveOperator &_A, RHS &_F,
                                                   bool _optimized_grow, int _assemble_matrix)
    : A(_A), F(_F), optimized_grow(_optimized_grow), assemble_matrix(_assemble_matrix),
      cA(A.cA), CA(A.CA), kappa(A.kappa),
      alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.),
      hms_galerkin(389), hms_residual(3079)
{
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;
    if (IsIndex1D<Index>::value) {
        //if (A.basis.d==2) {
        //    hms_galerkin = 6151;
        //    hms_residual = 24593;
        //}
        //else {
            hms_galerkin = 3079;
            hms_residual = 12289;
        //}
    }
    else if (IsIndex2D<Index>::value) {
        hms_galerkin = 3145739;
        hms_residual = 6291469;
    }
    else if (IsIndex3D<Index>::value) {
        hms_galerkin = 3145739;
        hms_residual = 6291469;
    }
    else {
        hms_galerkin = 3145739;
        hms_residual = 6291469;
    }
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
void
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::setParameters(T _alpha, T _omega, T _gamma, T _theta)
{
    alpha = _alpha;
    omega = _omega;
    gamma = _gamma;
    theta = _theta;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
Coefficients<Lexicographical,T,Index>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::SOLVE(T nuM1, T _eps, const char *filename,
                                               int NumOfIterations, T H1norm)
{
    T eps = _eps;
    Coefficients<Lexicographical,T,Index> w_k(hms_galerkin), g_kP1(hms_galerkin);
    IndexSet<Index>                       Lambda_kP1(hms_galerkin);
    T nu_kM1 = nuM1;
    T nu_k   = 0.;
    T total_time=0.;
    int numOfIterations=0;
    T timeApply=0.;
    T timeMatrixVector=0.;
    int maxIterations=100;
    int lengthOfResidual = 1;

    std::cerr << "GHS-ADWAV-SOLVE has started with the following parameters: " << std::endl;
    std::cerr << "  alpha=" << alpha << ", omega=" << omega << ", gamma=" << gamma << ", theta="
              << theta << std::endl;
    std::cerr << "  cA=" << cA << ", CA=" << CA << ", kappa=" << kappa << std::endl;

    std::ofstream file(filename);

    for (int i=1; i<=NumOfIterations; ++i) {
        Timer time;
        std::cerr << "*** " << i << ".iteration ***" << std::endl;
        time.start();
        std::cerr << "   GROW started." << std::endl;

        IndexSet<Index> Extension;
        Extension = this->GROW(w_k, theta*nu_kM1, nu_k, lengthOfResidual, timeApply);
        Lambda_kP1 += Extension;
        std::cerr << "   GROW finished." << std::endl;

        time.stop();
        total_time += time.elapsed();

        if (nu_k <=eps) break;

        Coefficients<Lexicographical,T,Index> Au, f, rhs;
        f = F(supp(w_k));
        T fu = w_k*f;
        Au = A.mv(supp(w_k),w_k);
        T uAu = w_k*Au;
        std::cerr.precision(16);
        T Error_H_energy = sqrt(fabs(std::pow(H1norm,(T)2.)- 2*fu + uAu));
        file << w_k.size() << " " << numOfIterations << " " << total_time << " " <<  nu_k << " "
                         << Error_H_energy << " " << timeApply << " " << timeMatrixVector << " "
                         << T(lengthOfResidual)/T(w_k.size())  << std::endl;

        time.start();
        std::cerr << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;

        //Attention: update in rhs1d is set! Linear complexity still holds with better convergence!

        //g_kP1 = F(Lambda_kP1);                // update rhs vector
        //rhs = F(gamma*nu_k);
        //g_kP1 = P(rhs,Lambda_kP1);  // compute with restriction, otherwise GROW does not work
        g_kP1 = F(gamma*nu_k);
        P(Lambda_kP1, g_kP1);
        T res;
        numOfIterations=CG_Solve(Lambda_kP1, A, w_k, g_kP1, res, gamma*nu_k,
                                 maxIterations,timeMatrixVector,assemble_matrix);
        std::cerr << "      ... required " << numOfIterations << " iterations." << std::endl;
        if (numOfIterations>maxIterations) {
            std::cerr << "   Attention: maximum number of iterations exceeded!" << std::endl;
        }
        //this->GALSOLVE(Lambda_kP1, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);
        std::cerr << "   GALSOLVE finished." << std::endl;

        nu_kM1 = nu_k;
        time.stop();
        total_time += time.elapsed();

        //std::stringstream coeff_filename;
        //coeff_filename << "adwav_coeff_" << w_k.size();
        //plotScatterCoeff2D(w_k, A.basis.first, A.basis.second, coeff_filename.str().c_str());

        std::cerr << std::endl;
    }
    return w_k;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
IndexSet<Index>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GROW(const Coefficients<Lexicographical,T,Index> &w,
                                              T nu_bar, T &nu, int &lengthOfResidual, T &timeApply)
{
    T zeta = 2.*(omega*nu_bar)/(1-omega);
    T r_norm = 0.;
    Coefficients<Lexicographical,T,Index> r(hms_residual);
    while (1) {
        r.setToZero();
        zeta *= 0.5;
        std::cerr << "      zeta = " << zeta << std::endl;
        Timer time_apply;
        time_apply.start();
        A.apply(w, 0.5*zeta, r);
        time_apply.stop();
        timeApply = time_apply.elapsed();
        std::cerr << "      Required time for APPLY: " << timeApply << std::endl;
        r -= F(0.5*zeta);

        r_norm = r.norm(2.);
        nu = r_norm + zeta;
        std::cerr << "      zeta = " << zeta << ", r_norm = " << r_norm << ", omega*r_norm = "
                  << omega*r_norm << ", nu = " << nu << std::endl;
        if (nu <= eps) break;
        if (zeta<=omega*r_norm) break;
    }
    lengthOfResidual = r.size();

    IndexSet<Index> Lambda;
    long double P_Lambda_r_norm_square = 0.0L;
    if (w.size() > 0) {
        for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
            //Lambda.insert((*it).first);
            P_Lambda_r_norm_square += std::pow(r[(*it).first],(T)2.);
            r.erase((*it).first);
        }
    }

    T threshbound = std::sqrt(1-alpha*alpha) * r.norm(2.)/std::sqrt(T(r.size()));
    Coefficients<Bucket,T,Index> r_bucket;
    std::cerr << "      size of r = " << r.size() << std::endl;

    r_bucket.bucketsort(r, threshbound);


    std::cerr << "      ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << alpha*r_norm << std::endl;
    if (nu > eps) {
        for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
            P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
            std::cerr << "      ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
            int addDoF = r_bucket.addBucketToIndexSet(Lambda,i);
            std::cerr << "      Added " << addDoF << " indices with norm " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
            if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) {
                if (optimized_grow) {
                    int addDoF = r_bucket.addBucketToIndexSet(Lambda,i+1);
                    std::cerr << "      Added " << addDoF << " indices with norm " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
                }
                break;
            }
        }
    }
    std::cerr << "      #Lambda_{k+1} / #Lambda_k = " << T(Lambda.size()+w.size())/T(w.size()) << std::endl;
    return Lambda;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
void
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GALSOLVE(const IndexSet<Index> &Lambda,
                                                  Coefficients<Lexicographical,T,Index> &w,
                                                  const Coefficients<Lexicographical,T,Index> &g,
                                                  T delta, T tol)
{
    if (assemble_matrix==1 || assemble_matrix==2) {
        std::map<Index,int,lt<Lexicographical,Index> >      row_indices;
        int d=A.basis.d;
        if (Lambda.size()==0) return;

        //Determine compression level
        int J=0;        //compression
        if      (d==2) {   J = -std::ceil(2*std::log(tol/((3*tol+3*delta)*kappa))); }
        else if (d==3) {   J = -std::ceil((1./1.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
        else if (d==4) {   J = -std::ceil((1./2.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
        else              { assert(0); }
        std::cerr << "      Estimated compression level for delta=" << delta << " and target tol=" << tol
                  << " : " << J << std::endl;

        //Assemble sparse matrix B
        unsigned long N = Lambda.size();
        //std::cerr << "    Assembling of B started with N=" << N << std::endl;


        flens::SparseGeMatrix<CRS<T,CRS_General> > B(N,N);
        Timer time_assemble;
        time_assemble.start();
        A.toFlensSparseMatrix(Lambda,Lambda,B,J);
        time_assemble.stop();
        std::cerr << "      Required time for matrix assembly: " << time_assemble.elapsed() << std::endl;

        Coefficients<Lexicographical,T,Index> r0, APPLY_Aw;
        Timer time_apply;
        time_apply.start();
        A.apply(w, tol/3.,Lambda, APPLY_Aw);
        time_apply.stop();
        std::cerr << "      Required time for APPLY: " << time_apply.elapsed() << std::endl;
        r0 = g - APPLY_Aw;

        DenseVector<Array<T> > rhs(N), x(N), res(N), Bx(N);

        const_coeff_it r0_end = r0.end();
        int row_count=1;
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
            const_coeff_it it = r0.find(*row);
            if (it != r0_end) rhs(row_count) = (*it).second;
            else              rhs(row_count) = 0.;
            ++row_count;
        }

        int iters = lawa::cg(B,x,rhs,tol/3.);
        Bx = B*x;
        res= Bx-rhs;
        T lin_res = std::sqrt(res*res);
        std::cerr << "      cg-method needed " << iters << " iterations, res=" << lin_res << std::endl;
        assert(lin_res<tol);


        const_coeff_it w_end = w.end();

        row_count=1;
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row) {

            w[*row] += x(row_count);
            ++row_count;
        }
    }

    else {

        T alpha, beta, rNormSquare, rNormSquarePrev;
        Coefficients<Lexicographical,T,Index> r, p;
        A.apply(w, tol/3.,Lambda, r);
        r -= g;
        p -= r;
        rNormSquare = r*r;
        for (int k=1; k<=1000; k++) {
            Coefficients<Lexicographical,T,Index> Ap;
            A.apply(p,tol/3.,Lambda, Ap);
            T pAp = p * Ap;
            //T pAp = A.innerproduct(p);
            alpha = rNormSquare/pAp;
            w += alpha*p;
            r += alpha*Ap;

            if (sqrt(rNormSquare)<=tol) {
                std::cerr << "      cg iterations: " << k << std::endl;
                break;
            }

            rNormSquarePrev = rNormSquare;
            rNormSquare = r*r;
            beta = rNormSquare/rNormSquarePrev;
            p *= beta;
            p -= r;
        }

    }
    return;
}

/*
        std::stringstream coefffile;
        coefffile << "adwav_coeffs_" << w_k.size();
        plotScatterCoeff2D(w_k, A.basis.first, A.basis.second, coefffile.str().c_str());
        */



/*
int count=0;
int sizeExtensionOfLambda=1;
for (const_coeff_abs_it it=r_abs.begin(); it!=r_abs.end(); ++it) {

    Lambda.insert((*it).second);
    ++count;

    P_Lambda_r_norm_square += std::pow((*it).first,2);
    std::cerr << "    Added " << (*it).second << ", " << (*it).first << ", now: ||P_{Lambda}r ||_2 = "
              << std::sqrt(P_Lambda_r_norm_square) << ", alpha*r_norm = "
              << alpha*r_norm << std::endl;

    if (count>=5*sizeExtensionOfLambda)                     break;
}
*/


        /*
        index_hashfunction<Index2D> hasher2d;
        std::map<size_t,Index2D> hash_values;
        int collisions = 0;
        T time_hash = 0.;
        for (const_coeff_it it=w_k.begin(); it!=w_k.end(); ++it) {
            time.start();
            size_t hash_value = hasher2d((*it).first);
            time.stop();
            time_hash += time.elapsed();
            int count = hash_values.count(hash_value);
            if (count==0) hash_values[hash_value] = (*it).first;
            else {
               std::cerr << "   Collision: " << (*it).first << " " << hash_value << " same as " << hash_values[hash_value] << " -> " << count << std::endl;
                ++collisions;
            }
        }
        std::cerr << "   Hash test: " << collisions << " / " << w_k.size() << " " << time_hash << std::endl;
        T percent_collisions = (T)collisions / (T) w_k.size();
        */
        /*
        std::stringstream coeff_filename;
        coeff_filename << "adwav_coeff_" << w_k.size();
        Coefficients<AbsoluteValue,T,Index1D> u_abs;
        u_abs = w_k;
        plotCoeff(u_abs, A.basis, coeff_filename.str().c_str());
        */
        /*
        std::stringstream coefffile;
        coefffile << "s_adwav_coeffs_" << i;
        plotScatterCoeff2D(w_k, A.basis.first, A.basis.second, coefffile.str().c_str());
        */

}   //namespace lawa
