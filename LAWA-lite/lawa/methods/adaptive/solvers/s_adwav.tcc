/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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
#include <lawa/aux/timer.h>
#include <limits.h>

namespace lawa {


template <typename T, typename Index, typename Basis, typename MA, typename RHS>
S_ADWAV<T,Index,Basis,MA,RHS>::S_ADWAV(const Basis &_basis, MA &_A, RHS &_F, T _contraction,
                                 T start_threshTol, T start_linTol, T start_resTol,
                                 int _NumOfIterations, int _MaxItsPerThreshTol, T _eps, int _MaxSizeLambda,
                                  T _resStopTol, std::vector<int> _Jmaxvec)
    : basis(_basis), A(_A), F(_F), contraction(_contraction), threshTol(start_threshTol), linTol(start_linTol),
      resTol(start_resTol), NumOfIterations(_NumOfIterations), MaxItsPerThreshTol(_MaxItsPerThreshTol), eps(_eps),
      MaxSizeLambda(_MaxSizeLambda), resStopTol(_resStopTol), Jmaxvec(_Jmaxvec)
{
    if(Jmaxvec.size() == 0){
        Jmaxvec.push_back(INT_MAX);
        Jmaxvec.push_back(INT_MAX);
        Jmaxvec.push_back(INT_MAX);
    }
    
    solutions.resize(NumOfIterations);
    residuals.resize(NumOfIterations);
    times.resize(NumOfIterations);
    toliters.resize(NumOfIterations);
    linsolve_iterations.resize(NumOfIterations);
}


template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve(const IndexSet<Index> &InitialLambda, const char *linsolvertype,
                                     const char *filename, int assemble_matrix, T H1norm)
{
    Timer timer;

    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable, DeltaLambda;
    Coefficients<Lexicographical,T, Index> u, f, Au, r;

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    T timeMatrixVector=0.;
    std::cout << "Simple adaptive solver started." << std::endl;

    std::ofstream file(filename);

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();

        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
        int iterations=0;

        if (strcmp(linsolvertype,"cg")==0) {
            iterations = CG_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol, 100, timeMatrixVector, assemble_matrix);
        }
        else if (strcmp(linsolvertype,"gmres")==0) {
            iterations = GMRES_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol,100, assemble_matrix);
        }
        else {
            assert(0);
            std::cerr << "Linsolver type " << linsolvertype << std::endl;
            exit(1);
        }
        linsolve_iterations[its] = iterations;
        std::cerr << "   ...finished with res=" << r_norm_LambdaActive << std::endl;

        // Attention: Relative threshold here removes lots of entries -> linear system was much
        // larger than necessary.
        T thresh_off_quot = 1./T(u.size());
        //u = THRESH(u,threshTol*u.norm(2.),false);
        //u = THRESH(u,threshTol,false, basis.d > 3 ? true : false);
        u = THRESH(u,threshTol,false, true);
        thresh_off_quot *= T(u.size());
        solutions[its] = u;
        LambdaThresh = supp(u);

        timer.stop();

        std::stringstream coeff_filename;
        coeff_filename << "s_adwav_coeff_" << u.size();
        //Coefficients<AbsoluteValue,T,Index1D> u_abs;
        //u_abs = u;
        //plotCoeff(u_abs, basis, coeff_filename.str().c_str());

        std::cout << "Computing error..." << std::endl;
        T time1 = timer.elapsed();

        Coefficients<Lexicographical,T,Index> Au_test;
        Au_test = A.mv(LambdaThresh,u);
        T uAu_test = u*Au_test;
        std::cerr << "    Estimated energy norm squared: " << uAu_test << std::endl;

        T Error_H_energy = 0.;
        if (H1norm>0) Error_H_energy = computeErrorInH1Norm(A, F, u, H1norm, true);
        std::cerr << "    Estimated energy error squared: " << Error_H_energy*Error_H_energy
                  << ", ref_value: " << H1norm*H1norm << std::endl;

        std::cout << "... finished." << std::endl;


        timer.start();
        //Computing residual
        std::cout << "   Computing DeltaLambda..." << std::endl;
        DeltaLambda = C(LambdaThresh, contraction, basis);

        std::cout << "... finished." << std::endl;
        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);



        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        //Au = mv(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            Au = mv_sparse(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;

        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;
        toliters[its] = linTol;

        //r = THRESH(r,threshTol*r.norm(2.),false,false);
        //r = THRESH(r,threshTol, false, basis.d > 3 ? true : false);
        r = THRESH(r,threshTol, false, true);

        LambdaActive = LambdaThresh+supp(r);


        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }

        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) {
            T total_time = time1+timer.elapsed();
            file << LambdaThresh.size() << " " << iterations << " " <<  total_time << " "
                 << estim_res << " " << Error_H_energy << " " << thresh_off_quot << std::endl;
            times[its] = total_time;
        }
        else {
            T total_time = times[its-1] + time1 + timer.elapsed();
            file << LambdaThresh.size() << " " << iterations << " "  << total_time << " "
                 << estim_res << " " << Error_H_energy << " " << thresh_off_quot << std::endl;
            times[its] = total_time;
        }

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cg-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;

        if(supp(u).size() > (unsigned int) MaxSizeLambda){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }
    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_cg(const IndexSet<Index> &InitialLambda, int assemble_matrix, T H1norm)
{
    Timer timer;

    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable, DeltaLambda;
    Coefficients<Lexicographical,T, Index> u, f, Au, r;

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    T timeMatrixVector=0.;

    std::cout << "Simple adaptive solver started." << std::endl;
    
    std::stringstream filename;
    filename << "s-adwav-otf.dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;
        
        timer.start();
        
        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        //std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = CG_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol, 1000, timeMatrixVector, assemble_matrix);
        linsolve_iterations[its] = iterations;
        //std::cout << "   ...finished." << std::endl;

        //Threshold step
        //std::cout << "Before THRESH: " << u << std::endl;
        //std::cout << "||u||_2 = " << u.norm(2.) << ", threshTol = " << threshTol << std::endl;
        // u = THRESH(u,threshTol*u.norm(2.),false);
        u = THRESH(u,threshTol,false, basis.d > 3 ? true : false);

				//int jmin,jmax;
        //getMinAndMaxLevel(LambdaActive, jmin, jmax);
        //std::cout << "Before THRESH: jmin = " << jmin << ", jmax = " << jmax << std::endl;
        //std::cout << "After THRESH: " << u << std::endl;

        solutions[its] = u;
        LambdaThresh = supp(u);
        //getMinAndMaxLevel(LambdaThresh, jmin, jmax);
        //std::cout << "After THRESH: jmin = " << jmin << ", jmax = " << jmax << std::endl;
        //std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;
        std::stringstream filename_coefficients;
        filename_coefficients << "coefficients_" << its+1;

        timer.stop();

        T time1 = timer.elapsed();
        T Error_H_energy = 0.;
        if (H1norm>0) Error_H_energy = computeErrorInH1Norm(A, F, u, H1norm);


        timer.start();
        //Computing residual
        DeltaLambda = C(LambdaThresh, contraction, basis, Jmaxvec[0], Jmaxvec[1]);            

        //std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        //std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        //std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        //Au = mv(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            Au = mv_sparse(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        //std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;
        toliters[its] = linTol;

        //r = THRESH(r,threshTol*r.norm(2.));
        r = THRESH(r,threshTol, false, basis.d > 3 ? true : false);
				LambdaActive = LambdaThresh+supp(r);


        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }

        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) {
            T total_time = time1+timer.elapsed();
            file << LambdaThresh.size() << " " << iterations << " " <<  total_time << " "
                 << estim_res << " " << Error_H_energy << std::endl;
            times[its] = total_time;
        }
        else {
            T total_time = times[its-1] + time1 + timer.elapsed();
            file << LambdaThresh.size() << " " << iterations << " "  << total_time << " "
                 << estim_res << " " << Error_H_energy << std::endl;
            times[its] = total_time;
        }

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cg-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;
        
        if((supp(u).size() > (unsigned int) MaxSizeLambda) || (estim_res < resStopTol)){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }
    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_cg_WO_XBSpline(const IndexSet<Index> &InitialLambda, int assemble_matrix, T H1norm)
{
    Timer timer;

    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable, DeltaLambda;
    Coefficients<Lexicographical,T, Index> u, f, Au, r;

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
		T timeMatrixVector = 0;
		
    std::cout << "Simple adaptive solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-realline-helmholtz-otf.dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();

        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = CG_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol, 1000, timeMatrixVector, assemble_matrix);
        std::cout << "   ...finished." << std::endl;

        //Threshold step
        //u = THRESH(u,threshTol);
        u = THRESH(u,threshTol, false, basis.d > 3 ? true : false);
				solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;

        timer.stop();
        T time1 = timer.elapsed();
        T Error_H_energy = 0.;
        if (H1norm>0) Error_H_energy = computeErrorInH1Norm(A, F, u, H1norm);

        timer.start();
        //Computing residual
        DeltaLambda = C_WO_XBSpline(LambdaThresh, contraction, basis);
        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        //Do not use mv_sparse before lambdaTilde works properly for very low levels!!
        //Au = mv_sparse(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            //Au = mv_sparse(DeltaLambda,A,u);
					Au = mv(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        file << LambdaThresh.size() << " " << estim_res << " " << Error_H_energy << std::endl;

        //r = THRESH(r,threshTol);
				r = THRESH(r,threshTol, false, basis.d > 3 ? true : false);
        LambdaActive = LambdaActive+supp(r);

        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) times[its] = time1+timer.elapsed();
        else        times[its] = times[its-1] + time1 + timer.elapsed();
        
        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cg-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;
                  
        if((supp(u).size() > (unsigned int) MaxSizeLambda) || (estim_res < resStopTol)){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }

    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_gmres(const IndexSet<Index> &InitialLambda, int assemble_matrix)
{
    Timer timer;
    
    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable, DeltaLambda;
    Coefficients<Lexicographical,T, Index> u, f, Au, r;

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-otf-gmres.dat";
    std::ofstream file(filename.str().c_str());
    T total_time = 0.;

    for (int its=0; its<NumOfIterations; ++its) {

    
        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);

        timer.start();

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   GMRES solver started with N = " << LambdaActive.size() << std::endl;
        int maxIterations = 1000;
        int iterations = GMRES_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol, maxIterations, assemble_matrix);
        std::cout << "   ...finished." << std::endl;



        //Threshold step
        // u = THRESH(u,threshTol);
        u = THRESH(u,threshTol, false, basis.d > 3 ? true : false);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;
        //int current_jmin, current_jmax;
        //getMinAndMaxLevel(LambdaThresh, current_jmin, current_jmax);
        //std::cout << "    Current minimal level: " << current_jmin << ", current maximal level: " << current_jmax << std::endl;

        //Computing residual
        DeltaLambda = C(LambdaThresh, contraction, basis, Jmaxvec[0], Jmaxvec[1]);            

        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;

        timer.stop();
        T time_galerkin = timer.elapsed();

        timer.start();
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;

        //Au = mv(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            Au = mv_sparse(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        //r = THRESH(r,threshTol);
        u = THRESH(u,threshTol, false, basis.d > 3 ? true : false);

        //LambdaActive = LambdaThresh+supp(r);
        LambdaActive = LambdaActive+supp(r);

        timer.stop();
        T time_res = timer.elapsed();

        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
        //if(old_res - estim_res < resTol){
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        times[its] = time_galerkin + time_res;

        total_time += times[its];
        file << LambdaThresh.size() << " " << total_time << " " << estim_res << std::endl;

        //times[its] = timer.elapsed();
        
        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", gmres-its = " << iterations;
        std::cout << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl;
        
        if((supp(u).size() > (unsigned int) MaxSizeLambda) || (estim_res < resStopTol)){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }
    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_gmresm(const IndexSet<Index> &InitialLambda, int assemble_matrix)
{
    Timer timer;
    
    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable, DeltaLambda;
    Coefficients<Lexicographical,T, Index> u, f, Au, r;
    
    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-otf-gmres.dat";
    std::ofstream file(filename.str().c_str());
    T total_time = 0.;
    
    for (int its=0; its<NumOfIterations; ++its) {
        
        
        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);
        
        timer.start();
        
        //Galerkin step
        T r_norm_LambdaActive = 0.0;
				int maxIterations = 1000;
        std::cout << "   GMRESM solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = GMRESM_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol, maxIterations, assemble_matrix);
        std::cout << "   ...finished." << std::endl;
        
        
        
        //Threshold step
        //u = THRESH(u,threshTol);
        u = THRESH(u,threshTol, false, basis.d > 3 ? true : false);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;
        //int current_jmin, current_jmax;
        //getMinAndMaxLevel(LambdaThresh, current_jmin, current_jmax);
        //std::cout << "    Current minimal level: " << current_jmin << ", current maximal level: " << current_jmax << std::endl;
        
        //Computing residual
        DeltaLambda = C(LambdaThresh, contraction, basis, Jmaxvec[0], Jmaxvec[1]);            

        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        
        timer.stop();
        T time_galerkin = timer.elapsed();
        
        timer.start();
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        
        //Au = mv(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            Au = mv_sparse(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;
        
        //r = THRESH(r,threshTol);
        r = THRESH(r,threshTol, false, basis.d > 3 ? true : false);
        //LambdaActive = LambdaThresh+supp(r);
        LambdaActive = LambdaActive+supp(r);
        
        timer.stop();
        T time_res = timer.elapsed();
        
        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            //if(old_res - estim_res < resTol){
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        times[its] = time_galerkin + time_res;
        
        total_time += times[its];
        file << LambdaThresh.size() << " " << total_time << " " << estim_res << std::endl;
        
        //times[its] = timer.elapsed();
        
        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", gmresm-its = " << iterations;
        std::cout << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl;
        
        if((supp(u).size() > (unsigned int) MaxSizeLambda) || (estim_res < resStopTol)){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }
    }
}


template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_cgls(const IndexSet<Index> &InitialLambda, int assemble_matrix)
{
    Timer timer;

    IndexSet<Index> LambdaActive, LambdaThresh, LambdaActivable;
    IndexSet<Index>      LambdaActive_test, DeltaLambda;
    Coefficients<Lexicographical,T, Index>      u, f, Au, r;

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive cgls time solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-cgls.dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
        std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();

        //Initialization step
        LambdaActive_test = LambdaActive + C_t(LambdaActive,contraction,basis);
        FillWithZeros(LambdaActive,u);
        f  = F (LambdaActive_test);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
				int maxIterations = 100000;
        int iterations = CGLS_Solve(LambdaActive_test, LambdaActive,
                                    A, u, f, r_norm_LambdaActive, linTol, maxIterations, assemble_matrix);
        std::cout << "   ...finished after " << iterations << " iterations with residual = " << r_norm_LambdaActive << std::endl;

        //Threshold step
        //u = THRESH(u,threshTol);
        u = THRESH(u,threshTol, false, basis.d > 3 ? true : false);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;

        //Computing residual for operator part
        DeltaLambda = LambdaThresh + C(LambdaThresh, contraction, basis);
        std::cout << "   Computing rhs for DeltaLambda_test (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda_test_operator..." << std::endl;
        //Au = mv_sparse(DeltaLambda,A,u);
        if (assemble_matrix==1) {
            Au = mv_sparse(DeltaLambda,A,u);
        }
        else {
            Au = A.mv(DeltaLambda,u);
        }
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        std::cout << "   ...finished" << std::endl;

        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        file << u.size() << " " << LambdaActive.size() << " " << estim_res << " " << iterations << " " << r_norm_LambdaActive << std::endl;

        //r = THRESH(r,threshTol);
        r = THRESH(r,threshTol, false, basis.d > 3 ? true : false);
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

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cgls-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;
                  
        if((supp(u).size() > (unsigned int) MaxSizeLambda) || (estim_res < resStopTol)){
            NumOfIterations = its+1;
            solutions.resize(NumOfIterations);
            residuals.resize(NumOfIterations);
            times.resize(NumOfIterations);
            toliters.resize(NumOfIterations);
            linsolve_iterations.resize(NumOfIterations);
            break;
        }

    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::set_parameters(T _contraction, T start_threshTol, T _linTol, 
                                              T _resTol, int _NumOfIterations, 
                                              int _MaxItsPerThreshTol, T _eps, int _MaxSizeLambda, 
                                              T _resStopTol, std::vector<int> _Jmaxvec)
{
    contraction = _contraction;
    threshTol = start_threshTol;
    linTol = _linTol;
    resTol = _resTol;
    NumOfIterations = _NumOfIterations;
    MaxItsPerThreshTol = _MaxItsPerThreshTol;
    eps = _eps;
    MaxSizeLambda = _MaxSizeLambda;
    resStopTol = _resStopTol;
    Jmaxvec = _Jmaxvec;
    if(Jmaxvec.size() == 0){
        Jmaxvec.push_back(INT_MAX);
        Jmaxvec.push_back(INT_MAX);
        Jmaxvec.push_back(INT_MAX);
    }
    
    solutions.resize(NumOfIterations);
    residuals.resize(NumOfIterations);
    times.resize(NumOfIterations);
    toliters.resize(NumOfIterations);
    linsolve_iterations.resize(NumOfIterations);
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::get_parameters(T& _contraction, T& _threshTol, T& _linTol, T& _resTol, 
                                              int& _NumOfIterations, int& _MaxItsPerThreshTol, T& _eps, 
                                              int& _MaxSizeLambda, T& _resStopTol, std::vector<int>& _Jmaxvec)
{
    _contraction = contraction;
    _threshTol = threshTol;
    _linTol = linTol;
    _resTol = resTol;
    _NumOfIterations = NumOfIterations;
    _MaxItsPerThreshTol = MaxItsPerThreshTol;
    _eps = eps;
    _MaxSizeLambda = MaxSizeLambda;
    _resStopTol = resStopTol;
    _Jmaxvec = Jmaxvec;
}





/*
std::stringstream coeff_filename;
coeff_filename << "s_adwav_coeff_" << u.size();
Coefficients<AbsoluteValue,T,Index1D> u_abs;
u_abs = u;
plotCoeff(u_abs, basis, coeff_filename.str().c_str());
*/
/*
std::stringstream coefffile;
coefffile << "s_adwav_coeffs_" << its;
plotScatterCoeff2D(u, A.basis.first, A.basis.second, coefffile.str().c_str());
*/


/*
       std::cout << "Computing eigenvalues..." << std::endl;
       T cB, CB;
       int N = LambdaActive.size();
       SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A_sparse(N,N);
       A.toFlensSparseMatrix(LambdaActive,LambdaActive,A_sparse);
       flens::DenseVector<flens::Array<T> > x(N);
       for (int i=1; i<=N; ++i) {
           x(i) = 1.;
       }
       std::cout << "powerMethod started." << std::endl;
       lawa::powerMethod(A_sparse,(T)1e-12,CB,x);
       std::cout << "powerMethod finished." << std::endl;
       for (int i=1; i<=N; ++i) {
           x(i) = 1.;
       }
       std::cout << "inversePowerMethod started." << std::endl;
       lawa::inversePowerMethod(A_sparse,(T)1e-12,cB,x);
       std::cout << "inversePowerMethod finished." << std::endl;

       std::cout << "  -> cB = " << cB << ", CB = " << CB << std::endl;
       */
}    //namespace lawa

