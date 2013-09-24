namespace lawa {

template <typename T, typename Basis, typename APPLY1D, typename RHS>
GHS_ADWAV1D<T,Basis,APPLY1D,RHS>::GHS_ADWAV1D(const Basis &_basis, APPLY1D &_Apply, RHS &_F)
    : basis(_basis), Apply(_Apply), A(Apply.A), F(_F),
      cA(0.), CA(0.), kappa(0.),
      alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.)
{
    cA    = Apply.parameters.cA;
    CA    = Apply.parameters.CA;
    kappa = Apply.parameters.kappa;
    Apply.parameters.getGHSADWAVParameters(alpha, omega, gamma, theta);


}

template <typename T, typename Basis, typename APPLY1D, typename RHS>
Coefficients<Lexicographical,T,Index1D>
GHS_ADWAV1D<T,Basis,APPLY1D,RHS>::SOLVE(T nuM1, T _eps, int NumOfIterations, T H1norm)
{
    T eps = _eps;
    Coefficients<Lexicographical,T,Index1D> w_k, w_kP1, g_kP1;
    IndexSet<Index1D> Lambda_kP1;
    T nu_kM1 = nuM1;
    T nu_k   = 0.;
    T total_time=0.;

    std::cerr << "GHS-ADWAV-SOLVE has started with the following parameters: " << std::endl;
    std::cerr << "  alpha=" << alpha << ", omega=" << omega << ", gamma=" << gamma << ", theta="
              << theta << std::endl;
    std::cerr << "  cA=" << cA << ", CA=" << CA << ", kappa=" << kappa << std::endl;

    std::stringstream filename;
    filename << "ghs-adwav-otf.dat";
    std::ofstream file(filename.str().c_str());

    for (int i=1; i<=NumOfIterations; ++i) {
        Timer time;
        std::cerr << "*** " << i << ".iteration ***" << std::endl;
        time.start();
        std::cerr << "  GROW started." << std::endl;
        Lambda_kP1 = GROW(w_k, theta*nu_kM1, nu_k);
        std::cerr << "  GROW finished." << std::endl;
        solutions.push_back(w_kP1);
        residuals.push_back(nu_k);
        times.push_back(total_time);
        if (nu_k <=eps) break;

        time.stop();
        total_time += time.elapsed();

        int current_jmin, current_jmax;
        getMinAndMaxLevel(supp(w_k), current_jmin, current_jmax);
        std::cout << "    Current minimal level: " << current_jmin << ", current maximal level: " << current_jmax << std::endl;
        T Error_H_energy = computeErrorInH1Norm(Apply.A, F, w_k, H1norm);
        file << w_k.size() << " " << total_time << " " <<  nu_k << " "
                         << Error_H_energy << std::endl;

        time.start();



        std::cerr << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;
        //g_kP1 = P(F(gamma*nu_k),Lambda_kP1);
        g_kP1 = F(Lambda_kP1);
        w_kP1 = GALSOLVE(Lambda_kP1, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);
        //T r_norm_LambdaActive;
        //int iterations = CG_Solve(Lambda_kP1, Apply.A, w_kP1, g_kP1, r_norm_LambdaActive, 1e-16);
        //std::cerr << "   iterations = " << iterations << ", residual = "
        //          << r_norm_LambdaActive << ", w_k = " << w_k << std::endl;
        std::cerr << "  GALSOLVE finished." << std::endl;

        nu_kM1 = nu_k;
        w_k = w_kP1;
        time.stop();
        total_time += time.elapsed();
        std::cerr << std::endl;
    }
    return w_k;
}

template <typename T, typename Basis, typename APPLY1D, typename RHS>
IndexSet<Index1D>
GHS_ADWAV1D<T,Basis,APPLY1D,RHS>::GROW(const Coefficients<Lexicographical,T,Index1D> &w,
                                       T nu_bar, T &nu)
{
    T zeta = 2.*(omega*nu_bar)/(1-omega);
    T r_norm = 0.;
    Coefficients<Lexicographical,T,Index1D> r, Aw, rhs;
    while (1) {
        zeta *= 0.5;
        rhs = F(0.5*zeta);
        Aw  = Apply(w, 0.5*zeta);
        r = rhs - Aw;
        r_norm = r.norm(2.);
        nu = r_norm + zeta;
        //std::cerr << "    zeta = " << zeta << ", r_norm = " << r_norm
        //          << ", omega*r_norm = " << omega*r_norm << ", nu = " << nu << std::endl;
        if (nu <= eps) break;
        if (zeta<=omega*r_norm) break;
    }


    Coefficients<AbsoluteValue,T,Index1D> r_abs;
    r_abs = r;
    IndexSet<Index1D> Lambda;
    T P_Lambda_r_norm_square = 0.;
    if (w.size() > 0) {
        Lambda = supp(w);
        for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
            P_Lambda_r_norm_square += std::pow(r[(*it).first],2.);
        }
    }

    //std::cerr << "   Before extension: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
    //          << ", alpha*r_norm = " << alpha*r_norm << std::endl;
    if (nu > eps) {
        for (const_coeff_abs_it it=r_abs.begin(); it!=r_abs.end(); ++it) {
            if (Lambda.count((*it).second) == 0) {
                Lambda.insert((*it).second);
                P_Lambda_r_norm_square += std::pow((*it).first,2);
                //std::cerr << "    Added " << (*it).second << ", now: ||P_{Lambda}r ||_2 = "
                //          << std::sqrt(P_Lambda_r_norm_square) << ", alpha*r_norm = "
                //          << alpha*r_norm << std::endl;
                if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) break;
            }
        }
    }
    return Lambda;
}

template <typename T, typename Basis, typename APPLY1D, typename RHS>
Coefficients<Lexicographical,T,Index1D>
GHS_ADWAV1D<T,Basis,APPLY1D,RHS>::GALSOLVE(const IndexSet<Index1D> &Lambda,
                                           const Coefficients<Lexicographical,T,Index1D> &g,
                                           const Coefficients<Lexicographical,T,Index1D> &w,
                                           T delta, T tol)
{
    int d=basis.d;
    Coefficients<Lexicographical,T,Index1D> ret;

    if (Lambda.size()==0) return ret;

    //Determine compression level
    int J=0;        //compression
    if         (d==2) {   J = -std::ceil(2*std::log(tol/((3*tol+3*delta)*kappa))); }
    else if (d==3) {   J = -std::ceil((1./1.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
    else              { assert(0); }
    //std::cerr << "   Estimated compression level for delta=" << delta << " and target tol=" << tol
    //          << " : " << J << std::endl;

    //Assemble sparse matrix B
    unsigned long N = Lambda.size();
    //std::cerr << "    Assembling of B started with N=" << N << std::endl;
    flens::SparseGeMatrix<CRS<T,CRS_General> > B(N,N);
    std::map<Index1D,int,lt<Lexicographical,Index1D> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    Apply.A.compression.setParameters(Lambda);
    for (const_set_it col=Lambda.begin(); col!=Lambda.end(); ++col, ++col_count) {
        IndexSet<Index1D> LambdaRowSparse = Apply.A.compression.SparsityPattern(*col, Lambda, J);
        for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
            T tmp = Apply.A(*row,*col);
            if (fabs(tmp)>0)    B(row_indices[*row],col_count) = 0.5*tmp;
            if (fabs(tmp)>0)    B(col_count,row_indices[*row]) = 0.5*tmp;
        }
    }
    B.finalize();
    //std::cerr << "    Assembling of B finished." << std::endl;

    //Compute approximate initial residual
    //std::cerr << "    Solving Bx=r0." << std::endl;
    Coefficients<Lexicographical,T,Index1D> r0, APPLY_Aw;
    APPLY_Aw = Apply(w, tol/3.);
    r0 = g - P(APPLY_Aw, Lambda);

    DenseVector<Array<T> > rhs(N), x(N), res(N), Bx(N);
    row_count=1;
    const_coeff_it r0_end = r0.end();
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
        const_coeff_it it = r0.find(*row);
        if (it != r0_end) rhs(row_count) = (*it).second;
        else                      rhs(row_count) = 0.;
    }
    //std::cerr << "    cg-method started with rhs " << rhs << std::endl;
    int iters = lawa::cg(B,x,rhs,tol);
    linsolve_iterations.push_back(iters);
    Bx = B*x;
    res= Bx-rhs;
    T lin_res = std::sqrt(res*res);
    std::cerr << "    cg-method needed " << iters << " iterations, res=" << lin_res << std::endl;
    assert(lin_res<tol);
    row_count = 1;
    const_coeff_it w_end = w.end();
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
        T tmp = 0.;
        const_coeff_it it = w.find(*row);
        if (it != w_end) tmp = (*it).second;
        ret[*row] = tmp+x(row_count);
    }
    //std::cerr << "    Calculated Bx=r0, w+x=" << ret << std::endl;
    return ret;
}


}    //namespace lawa

