#include <lawa/aux/timer.h>

namespace lawa {


template <typename T, typename Index, typename SpaceIndex, typename Basis, typename MA, typename RHSOperator, typename RHSInitialCond>
S_ADWAV_SPACETIME<T,Index,SpaceIndex,Basis,MA,RHSOperator,RHSInitialCond>::S_ADWAV_SPACETIME
                                                     (const Basis &_basis, MA &_A, RHSOperator &_F, RHSInitialCond &_U0,
                                                      T _contraction, T _threshTol, T _linTol, T _resTol,
                                                      int _NumOfIterations, int _MaxItsPerThreshTol, T _eps)
    : basis(_basis), A(_A), F(_F), U0(_U0), contraction(_contraction), threshTol(_threshTol), linTol(_linTol),
      resTol(_resTol), NumOfIterations(_NumOfIterations), MaxItsPerThreshTol(_MaxItsPerThreshTol), eps(_eps)
{
    solutions.resize(NumOfIterations);
    residuals.resize(NumOfIterations);
    times.resize(NumOfIterations);
}


template <typename T, typename Index, typename SpaceIndex, typename Basis, typename MA, typename RHSOperator, typename RHSInitialCond>
void
S_ADWAV_SPACETIME<T,Index,SpaceIndex,Basis,MA,RHSOperator,RHSInitialCond>::solve_cgls(const IndexSet<Index> &InitialLambda)
{
    Timer timer;

    int d=InitialLambda.d, d_=InitialLambda.d_;
    IndexSet<Index> LambdaActive(d,d_), LambdaThresh(d,d_), LambdaActivable(d,d_);
    IndexSet<Index>      LambdaActive_test_operator(d,d_), DeltaLambda_test_operator(d,d_);
    IndexSet<SpaceIndex> LambdaActive_test_initcond(d,d_), DeltaLambda_test_initcond(d,d_);
    Coefficients<Lexicographical,T, Index>      u(d,d_), f(d,d_), Au(d,d_), r(d,d_);
    Coefficients<Lexicographical,T, SpaceIndex> u0(d,d_), r0(d,d_), Au_0(d,d_);

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive space time solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-st-realline-heat1d_" << d << "_" << d_ << ".dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();


        //Initialization step
        LambdaActive_test_operator = LambdaActive + C_t(LambdaActive,contraction,basis);
        LambdaActive_test_initcond = extractSpaceIndices(LambdaActive);
        std::cout << "   Number of space indices: " << LambdaActive_test_initcond.size() << std::endl;
        FillWithZeros(LambdaActive,u);
        f  = F (LambdaActive_test_operator);
        u0 = U0(LambdaActive_test_initcond);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CGLS solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = CGLS_Solve(LambdaActive_test_operator, LambdaActive_test_initcond, A, LambdaActive,
                                    u, f, u0, r_norm_LambdaActive, linTol, 100000);
        std::cout << "   ...finished after " << iterations << " iterations with residual = " << r_norm_LambdaActive << std::endl;

        //Threshold step
        u = THRESH(u,threshTol);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;

        //Computing residual for operator part
        DeltaLambda_test_operator = LambdaThresh + C(LambdaThresh, contraction, basis);
        std::cout << "   Computing rhs for DeltaLambda_test_operator (size = " << DeltaLambda_test_operator.size() << ")" << std::endl;
        f = F(DeltaLambda_test_operator);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda_test_operator = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda_test_operator..." << std::endl;
        Au = mv_sparse(DeltaLambda_test_operator,A,u);
        r  = Au-f;
        T r_norm_DeltaLambda_test_operator = r.norm(2.);
        std::cout << "   ...finished" << std::endl;

        DeltaLambda_test_initcond = extractSpaceIndices(DeltaLambda_test_operator);
        u0 = U0(DeltaLambda_test_initcond);
        T u0_norm_DeltaLambda_test_initcond = u0.norm(2.);
        Au_0= mv(DeltaLambda_test_initcond,A,u);
        r0 = Au_0 - u0;
        T r0_norm_DeltaLambda_initcond = r0.norm(2.);

        T numerator   = r_norm_DeltaLambda_test_operator*r_norm_DeltaLambda_test_operator +
                        r0_norm_DeltaLambda_initcond*r0_norm_DeltaLambda_initcond;
        T denominator = f_norm_DeltaLambda_test_operator*f_norm_DeltaLambda_test_operator +
                        u0_norm_DeltaLambda_test_initcond*u0_norm_DeltaLambda_test_initcond;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        file << u.size() << " " << LambdaActive.size() << " " << estim_res << " " << iterations << " " << r_norm_LambdaActive << std::endl;

        r = THRESH(r,threshTol);
        LambdaActive = LambdaThresh+supp(r);

        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            //resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) times[its] = timer.elapsed();
        else        times[its] = times[its-1] + timer.elapsed();

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cg-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;

    }
}

}    //namespace lawa

