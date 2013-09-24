namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GHS_NONSYM_ADWAV(AdaptiveOperator &_A, RHS &_F, PP_AdaptiveOperator &_PP_A, PP_RHS &_PP_F,
                   bool _optimized_grow, bool _assemble_matrix)
    : A(_A), F(_F), PP_A(_PP_A), PP_F(_PP_F),
      optimized_grow(_optimized_grow), assemble_matrix(_assemble_matrix),
      cA(A.cA), CA(A.CA), kappa(A.kappa),
      alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.)
{
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
void
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::setParameters(T _alpha, T _omega, T _gamma, T _theta)
{
    alpha = _alpha;
    omega = _omega;
    gamma = _gamma;
    theta = _theta;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
Coefficients<Lexicographical,T,Index>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::SOLVE(T nuM1, T _eps, const char *filename, int NumOfIterations, T H1norm)
{
    T eps = _eps;
    Coefficients<Lexicographical,T,Index> w_k(SIZEHASHINDEX2D), g_kP1(SIZEHASHINDEX2D);
    IndexSet<Index> Lambda_kP1;
    T nu_kM1 = nuM1;
    T nu_k   = 0.;
    T total_time=0.;
    T timeApply=0.;
    int numOfIterations=0;
    T timeMatrixVector=0.;
    int maxIterations=300;
    int lengthOfResidual = 1;


    std::cerr << "GHS-NONSYM-ADWAV-SOLVE has started with the following parameters: " << std::endl;
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
        Lambda_kP1 = Lambda_kP1 + Extension;

        std::cerr << "   GROW finished." << std::endl;

        time.stop();
        total_time += time.elapsed();

        if (nu_k <=eps) break;

        /*
        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(Lambda_kP1, Lambda_x, Lambda_y);
        int jmin_x, jmax_x, jmin_y, jmax_y;
        getMinAndMaxLevel(Lambda_x, jmin_x, jmax_x);
        getMinAndMaxLevel(Lambda_y, jmin_y, jmax_y);
        std::cerr << "   Current jmax  = (" << jmax_x << ", " << jmax_y << ")" << std::endl;
        */
        Coefficients<Lexicographical,T,Index> PP_Au, PP_f, u;
        u = w_k;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            u[(*it).first] *= A.prec((*it).first);
            u[(*it).first] *= 1./PP_A.prec((*it).first);
        }
        //Maybe, we need this...
        //PP_f = F(supp(u));
        PP_f = PP_F(supp(u));
        T fu = u*PP_f;
        PP_Au = PP_A.mv(supp(u), u);
        T uAu = u*PP_Au;
        T Error_H_energy = sqrt(fabs(std::pow(H1norm,2.)- 2*fu + uAu));
        Coefficients<Lexicographical,T,Index> Au_M_f;
        A.apply(w_k,0.,Au_M_f,NoTrans);
        Au_M_f -= F(0.0001*nu_k);
        file << w_k.size() << " " << numOfIterations << " " << total_time << " " <<  nu_k << " "
                         << Error_H_energy << " " << Au_M_f.norm(2.) << " "  << T(lengthOfResidual)/T(w_k.size())  << std::endl;



        time.start();
        std::cerr << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;

        /*
         * Attention: update in rhs1d is set! Linear complexity still holds with better convergence!
         */
        Coefficients<Lexicographical,T,Index> rhs;
        g_kP1 = F(Lambda_kP1);                // update rhs vector
        std::cerr << "      Tolerance for GALSOLVE rhs: " << gamma*nu_k*(1./A.CA) << std::endl;
        g_kP1 = F(gamma*nu_k*(1./A.CA));
        //g_kP1 = P(rhs,Lambda_kP1);
        // compute with restriction, otherwise GROW does not work
        // "peaks" in the error estimator may arise from omega too small: if the approximation of rhs
        // for solving the Galerkin system is much smaller than the corresponding tolerance in GROW,
        // such effects can show up.
        T res;
        numOfIterations=CGLS_Solve(Lambda_kP1, A, w_k, g_kP1, res, gamma*nu_k,
                                   maxIterations,assemble_matrix);
        std::cerr << "      ... required " << numOfIterations << " iterations." << std::endl;
        //this->GALSOLVE(Lambda_kP1, Extension, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);
        std::cerr << "   GALSOLVE finished." << std::endl;

        nu_kM1 = nu_k;
        time.stop();
        total_time += time.elapsed();
        std::cerr << std::endl;
    }
    return w_k;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
IndexSet<Index>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GROW(const Coefficients<Lexicographical,T,Index> &w, T nu_bar, T &nu, int &lengthOfResidual,
       T &timeApply)
{
    T zeta = 2.*(omega*nu_bar)/(1-omega);
    T r_norm = 0.;
    Coefficients<Lexicographical,T,Index> r(SIZEHASHINDEX2D), AtAw(SIZEHASHINDEX2D),
                                          help(SIZEHASHINDEX2D);

    A.apply(w, zeta, help, cxxblas::NoTrans);
    A.apply(help, zeta, AtAw, cxxblas::Trans);
    std::cerr << "      size of A^T A w: " << AtAw.size() << std::endl;
    help.setToZero();
    while (1) {
        zeta *= 0.5;

        std::cerr << "      zeta = " << zeta << std::endl;
        help = F(zeta);
        A.apply(help, zeta, r, cxxblas::Trans);
        std::cerr << "      size of A^T f:    " << r.size() << std::endl;
        r -= AtAw;

        r_norm = r.norm(2.);
        nu = r_norm + zeta*A.CA;
        std::cerr << "      zeta = " << zeta << ", r_norm = " << r_norm << ", omega*r_norm = "
                  << omega*r_norm << ", nu = " << nu << std::endl;
        if (nu <= eps) break;
        if (zeta*A.CA<=omega*r_norm) break;
    }
    lengthOfResidual = r.size();

    IndexSet<Index> Lambda;
    long double P_Lambda_r_norm_square = 0.0L;
    if (w.size() > 0) {
        for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
            //Lambda.insert((*it).first); -> Extension is calculated!!
            P_Lambda_r_norm_square += std::pow(r[(*it).first],2.);
            r.erase((*it).first);
        }
    }

    T threshbound = std::sqrt(1-alpha*alpha) * r.norm(2.)/std::sqrt(T(r.size()));
    Coefficients<Bucket,T,Index> r_bucket;


    r_bucket.bucketsort(r, threshbound);
    //std::cerr << r_bucket << std::endl;

    std::cout << "      size of r = " << r.size() << std::endl;
    std::cerr << "      Before extension: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << alpha*r_norm << std::endl;
    if (nu > eps) {
        int addDoF = 0;
        for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
            P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
            addDoF += r_bucket.addBucketToIndexSet(Lambda,i);
            std::cerr << "      Currently " << addDoF << " added indices, ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
            if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) {
                if (optimized_grow) {
                    //while ((T) addDoF /(T) w.size() < 0.07) {
                        ++i;
                        P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
                        addDoF += r_bucket.addBucketToIndexSet(Lambda,i);
                        std::cerr << "      Currently " << addDoF << " added indices, ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
                    //}
                }
                break;
            }
        }
    }
    r.clear();
    AtAw.clear();
    help.clear();
    //getchar();
    return Lambda;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
void
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GALSOLVE(const IndexSet<Index> &Lambda, const IndexSet<Index> &Extension,
           const Coefficients<Lexicographical,T,Index> &g,
           Coefficients<Lexicographical,T,Index> &w, T /*delta*/, T tol)
{
    IndexSet<Index> LambdaCol, LambdaRow;
    LambdaCol = Lambda + Extension;
    Coefficients<Lexicographical,T,Index> help1, help2, help3;
    for (const_set_it it=LambdaCol.begin(); it!=LambdaCol.end(); ++it) {
        help1[*it] = 1.;
    }
    A.apply(help1,0.,help2,NoTrans);
    LambdaRow = supp(help2);
    if (assemble_matrix) {


        flens::SparseGeMatrix<CRS<T,CRS_General> > B(LambdaRow.size(),LambdaCol.size());
        A.toFlensSparseMatrix(LambdaRow,LambdaCol,B,1);

        std::cerr << "      LambdaRow.size() = " << LambdaRow.size() << ", LambdaCol.size() = " << LambdaCol.size() << std::endl;

        DenseVector<Array<T> > g_vec(LambdaRow.size()), x_vec(LambdaCol.size()), Ax_vec(LambdaRow.size());
        const_coeff_it g_end = g.end();
        const_coeff_it w_end = w.end();
        int row_count=1;
        for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row) {
            const_coeff_it it = g.find(*row);
            if (it != g_end) g_vec(row_count) = (*it).second;
            else             g_vec(row_count) = 0.;
            ++row_count;
        }
        row_count=1;
        for (const_set_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row) {
            const_coeff_it it2 = w.find(*row);
            if (it2 != w_end) x_vec(row_count) = (*it2).second;
            else              x_vec(row_count) = 0.;
            ++row_count;
        }
        std::cerr << "      target tolerance: " << tol << std::endl;
        int iters = lawa::cgls(B,x_vec,g_vec,tol);

        std::cerr << "      number of iterations: " << iters << ", res = " << tol << std::endl;

        Coefficients<Lexicographical,T,Index> x;
        row_count=1;
        for (const_set_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row) {
            w[(*row)] = x_vec(row_count);
            ++row_count;
        }
    }
    else {
        T alpha_cgls, beta_cgls, gammaPrev_cgls, gamma_cgls;
        Coefficients<Lexicographical,T,Index> r, q, s, p;

        A.apply(w,0.,r,NoTrans);
        r -= g;
        r *= -1.;
        A.apply(r,0.,LambdaCol,s,Trans);
        p = s;
        gammaPrev_cgls = s*s;
        std::cerr << "      gammaPrev = " << gammaPrev_cgls << std::endl;

        for (int k=1; k<=1000; k++) {
           q.setToZero();
           A.apply(p,0.,q,NoTrans);   //q = A*p;
           alpha_cgls = gammaPrev_cgls/(q*q);
           w +=   alpha_cgls*p;
           r -=   alpha_cgls*q;
           s.setToZero();
           A.apply(r,0.,LambdaCol,s,Trans);  // flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);

           gamma_cgls = s*s;

           std::cerr << "     iteration: " << k << " : current error = " << sqrt(gamma_cgls) << " (tol=" << tol << ")" << std::endl;
           if (sqrt(gamma_cgls)<=tol) {
               std::cerr << "      inner iterations: " << k << std::endl;
               break;
           }

           beta_cgls  = gamma_cgls/gammaPrev_cgls;
           p *= beta_cgls;
           p += s;
           gammaPrev_cgls = gamma_cgls;
        }

    }
}

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

}   // namespace lawa
