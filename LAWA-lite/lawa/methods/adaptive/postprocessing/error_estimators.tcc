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
 #include <iomanip>
 #include <lawa/righthandsides/righthandsides.h>
 #include <lawa/operators/operators.h>
 
namespace lawa {

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_Au_M_f(MA &A, RHS &F, const Coefficients<Lexicographical,T,Index> & u,
                     const IndexSet<Index> &LambdaRow)
{
    Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f(u.d,u.d_), res(u.d,u.d_);
    Au = mv_sparse(LambdaRow,A,u);
    f  = F(LambdaRow);
    res = Au-f;
    return res.norm(2.)/f.norm(2.);
}

template <typename T, typename Index, typename MA, typename RHS>
T
computeErrorInH1Norm(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
                     T HNormOfExactSolution, bool optimized)
{
    Coefficients<Lexicographical,T,Index> Au, f_H;
    IndexSet<Index> LambdaRow;
    LambdaRow = supp(u);
    if (optimized) {
        Au = A_H.mv(LambdaRow,u);
    }
    else {
        Au = mv_sparse(LambdaRow,A_H,u);
    }
    //Au = mv(supp(u),A_H,u);
    //std::cerr << "Attention: Slow matrix vector product..." << std::endl;
    T uAu = u*Au;
    f_H   = F_H(supp(u));
    T fu  = f_H*u;
    std::cerr << "   Estimated Energy norm squared: " << uAu << std::endl;

    return std::sqrt(fabs(std::pow(HNormOfExactSolution,(T)2.)- 2*fu + uAu));

}

template <typename T, typename Index, typename SOLVER, typename MA_H, typename RHS_H>
void
postprocessing_H1(SOLVER& Solver, MA_H &A_H1, RHS_H &F_H1, T H1norm, const char* filename)
{

    Coefficients<Lexicographical,T,Index> u;
    std::ofstream file(filename);

    for (int i=0; i<int(Solver.solutions.size()); ++i) {
        u = Solver.solutions[i];
        T ErrorH1Norm = computeErrorInH1Norm(A_H1, F_H1, u, H1norm);
        file      << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
        std::cerr << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
    }
    file.close();
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_L2(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += prec * u[*it] * prec * u[*it];
        }
    }
    return std::sqrt(error_est);
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_H1(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u[*it] * prec * u[*it];
        }
    }
    return std::sqrt(error_est);
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_W0T(Coefficients<Lexicographical,T,Index2D> & u, 
                            Coefficients<Lexicographical,T,Index2D> & u_exact,
                            const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    
    for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
        T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        T prec = P((*it));
        if( ExpandedLambda.count(*it) == 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u[*it] * prec * u[*it];
        }
    }
    
    return std::sqrt(error_est);
}

} // namespace lawa

