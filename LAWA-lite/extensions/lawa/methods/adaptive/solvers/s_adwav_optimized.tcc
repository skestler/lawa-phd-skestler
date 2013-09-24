namespace lawa {


template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
S_ADWAV_Optimized<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>::S_ADWAV_Optimized
(AdaptiveOperator &_A, RHS &_F, PP_AdaptiveOperator &_PP_A, PP_RHS &_PP_F,
 T _contraction, T start_threshTol, T start_linTol,
 T start_resTol, int _NumOfIterations, int _MaxItsPerThreshTol, T _eps)
: A(_A), F(_F), PP_A(_PP_A), PP_F(_PP_F),
  contraction(_contraction), threshTol(start_threshTol),
  linTol(start_linTol), resTol(start_resTol), NumOfIterations(_NumOfIterations),
  MaxItsPerThreshTol(_MaxItsPerThreshTol), eps(_eps)
{
    residuals.resize(NumOfIterations);
    times.resize(NumOfIterations);
    toliters.resize(NumOfIterations);
    linsolve_iterations.resize(NumOfIterations);
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
void
S_ADWAV_Optimized<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>::
solve(const IndexSet<Index> &InitialLambda, Coefficients<Lexicographical,T,Index> &u,
      const char *linsolvertype, const char *filename, int assemble_matrix, T H1norm)
{
    Timer timer;

    IndexSet<Index> Lambda(SIZEHASHINDEX2D), LambdaThresh(SIZEHASHINDEX2D), DeltaLambda(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T, Index> f(SIZEHASHINDEX2D), r(SIZEHASHINDEX2D);

    Lambda = InitialLambda;
    T old_res = linTol;
    int its_per_threshTol=0;
    T timeMatrixVector=0.;
    std::cout << "Simple adaptive solver started." << std::endl;

    std::ofstream file(filename);

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();

        //Initialization step
        FillWithZeros(Lambda,u);
        f = F(Lambda);
        T f_norm_Lambda = f.norm(2.);

        //Galerkin step
        T r_norm_Lambda = 0.0;
        std::cout << "   Solver " << linsolvertype << " started with N = " << Lambda.size() << std::endl;
        int iterations=0;

        if (strcmp(linsolvertype,"cg")==0) {
            iterations = CG_Solve(Lambda, A, u, f, r_norm_Lambda, linTol, 100, timeMatrixVector, assemble_matrix);
        }
        else if (strcmp(linsolvertype,"gmres")==0) {
            iterations = GMRES_Solve(Lambda, A, u, f, r_norm_Lambda, std::min(1e-2,0.5*old_res),35, assemble_matrix);
        }
        else if (strcmp(linsolvertype,"cgls")==0) {
            iterations = CGLS_Solve(Lambda, A, u, f, r_norm_Lambda, std::min(1e-2,0.5*old_res),1000, assemble_matrix);
        }
        else {
            assert(0);
            std::cerr << "Linsolver type " << linsolvertype << std::endl;
            exit(1);
        }
        linsolve_iterations[its] = iterations;
        std::cerr << "   ...finished with res=" << r_norm_Lambda << std::endl;
        // Attention: Relative threshold here removes lots of entries
        T thresh_off_quot = 1./T(u.size());


        //u = THRESH(u,0.5*threshTol,false, A.basis.d > 3 ? true : false);
        if (its!=NumOfIterations) {
            //u = THRESH(u,0.5*threshTol,false, A.basis.d > 3 ? true : false);
            u = THRESH(u,0.5*threshTol,false, true);
            Lambda = supp(u);
        }
        thresh_off_quot *= T(u.size());

        int N = Lambda.size();

        timer.stop();
        T time1 = timer.elapsed();

        std::cout << "Computing error..." << std::endl;
        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(Lambda, Lambda_x, Lambda_y);
        int jmin_x, jmax_x, jmin_y, jmax_y;
        getMinAndMaxLevel(Lambda_x, jmin_x, jmax_x);
        getMinAndMaxLevel(Lambda_y, jmin_y, jmax_y);
        std::cerr << "   Current jmax  = (" << jmax_x << ", " << jmax_y << ")" << std::endl;
        T Error_H_energy=0.;
        if (H1norm>0) {
            Coefficients<Lexicographical,T,Index> PP_Au, PP_f, u_tmp;
            u_tmp = u;
            for (const_coeff_it it=u_tmp.begin(); it!=u_tmp.end(); ++it) {
                u_tmp[(*it).first] *= A.prec((*it).first);
                u_tmp[(*it).first] *= 1./PP_A.prec((*it).first);
            }
            PP_f = PP_F(supp(u_tmp));
            T fu = u_tmp*PP_f;
            PP_Au = PP_A.mv(supp(u_tmp), u_tmp);
            T uAu = u_tmp*PP_Au;
            std::cerr << "   Estim. energy norm squared: " << u*PP_Au << ", ref. value: " << std::pow(H1norm,2.) << std::endl;
            Error_H_energy = sqrt(fabs(std::pow(H1norm,2.)- 2*fu + uAu));
        }
//        if (H1norm>0) Error_H_energy = computeErrorInH1Norm(A, F, u, H1norm, true);
        std::cout << "... finished." << std::endl;


        timer.start();
        //Computing residual
        r.setToZero();
        std::cout << "   Computing DeltaLambda..." << std::endl;
        DeltaLambda = C(Lambda, contraction, A.basis);
        std::cout << "... finished." << std::endl;
        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        //Au = mv(DeltaLambda,A,u);

        r.setToZero();
        A.apply(u, threshTol, DeltaLambda, r);
        r-=f;

        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = std::pow(r_norm_DeltaLambda,2.) + std::pow(r_norm_Lambda,2.);
        T denominator = std::pow(f_norm_DeltaLambda,2.) + std::pow(f_norm_Lambda,2.);
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;
        toliters[its] = linTol;

        //r = THRESH(r,threshTol*r.norm(2.),false,false);
        // f is currently not need and is overwritten in the next iteration anyway
        //f = THRESH(r, threshTol, false, A.basis.d > 3 ? true : false);
        f = THRESH(r, threshTol, false, true);

        Lambda +=  supp(f);

        timer.stop();
        if (its==0) {
            T total_time = time1+timer.elapsed();
            file << N << " " << iterations << " " <<  total_time << " "
                 << estim_res << " " << Error_H_energy << " " << thresh_off_quot << " " << threshTol << " " << linTol << std::endl;
            times[its] = total_time;
        }
        else {
            T total_time = times[its-1] + time1 + timer.elapsed();
            file << N << " " << iterations << " "  << total_time << " "
                 << estim_res << " " << Error_H_energy << " " << thresh_off_quot << " " << threshTol << " " << linTol << std::endl;
            times[its] = total_time;
        }

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size()
                  << ", cg-its = " << iterations << ", residual = " << estim_res
                  << " , current threshTol = " << threshTol << std::endl << std::endl;



        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }

        ++its_per_threshTol;
        old_res = estim_res;
    }
    f.clear();
    r.clear();
    return;
}

}   // namespace lawa
