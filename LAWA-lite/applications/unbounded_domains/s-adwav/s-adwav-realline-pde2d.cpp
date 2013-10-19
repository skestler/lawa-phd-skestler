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

#include <iostream>
#include <unistd.h>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>

typedef double T;
using namespace lawa;
using namespace std;

//Iterator definitions
typedef IndexSet<Index2D>::const_iterator                                 const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator           const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator             const_coeff2d_abs_it;

//Basis definitions
typedef Basis<T,Primal,R,SparseMulti>                                     SparseMW_Basis1D;
typedef TensorBasis2D<Adaptive, SparseMW_Basis1D,SparseMW_Basis1D>        SparseMW_Basis2D;

//Operator definitions
typedef AdaptivePDEOperatorOptimized2D<T,Primal,R,SparseMulti,
                                        Primal,R,SparseMulti>             SparseMW_MA;
//typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,R,SparseMulti,
//                                         Primal,R,SparseMulti>          SparseMW_MA;
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,R,SparseMulti,
                                             Primal,R,SparseMulti>        SparseMW_H1_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_MA>     SparseMW_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_H1_MA>  SparseMW_H1_Prec;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,SparseMW_Basis2D >                               SparseMW_SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SparseMW_SeparableRhsIntegral2D,
                             SparseMW_SeparableRhsIntegral2D>             SparseMW_SumOfSeparableRhsIntegral2D;

typedef RHS<T,Index2D,SparseMW_SumOfSeparableRhsIntegral2D,
            SparseMW_Prec>                                                SparseMW_SumOfSeparableRhs;
typedef RHS<T,Index2D,SparseMW_SumOfSeparableRhsIntegral2D,
            SparseMW_H1_Prec>                                             SparseMW_H1_SumOfSeparableRhs;

//Algorithm definition
typedef S_ADWAV_Optimized<T,Index2D,SparseMW_MA,SparseMW_SumOfSeparableRhs,
                          SparseMW_H1_MA,SparseMW_H1_SumOfSeparableRhs>   SparseMW_S_ADWAV_SOLVER_Optim_SeparableRhs;

int main (int argc, char *argv[]) {
    if (argc!=8) {
        cout << "usage " << argv[0] << " basistype d d_ jmin_x jmin_y example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);
    int example=atoi(argv[6]);
    int NumOfIterations=atoi(argv[7]);

    if (((JMINOFFSET + j0_x) < 0) || ((JMINOFFSET + j0_y) < 0)) {
        cout << "Please re-adjust #JMINOFFSET in index.h" << endl;
        exit(1);
    }

    T reaction=1., convection_x=0., convection_y=0., diffusion_y=1.;
    T contraction = 0.125;
    T threshTol = 0.4;
    T cgTol = 0.1*threshTol;//1e-12;
    T resTol=1e-4;

    IndexSet<Index2D> InitialLambda;
    Index1D index_x(j0_x,0,XBSpline);
    Index1D index_y(j0_y,0,XBSpline);
    InitialLambda.insert(Index2D(index_x,index_y));

    stringstream convfilename;
    convfilename << "s_adwav_conv_realline_pde2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_cx_" << convection_x << "_cy_"
                 << convection_y << "_ay_" << diffusion_y << "_" << argv[6] << ".dat";
    stringstream plotfilename;
    plotfilename << "s_adwav_plot_realline_pde2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_cx_" << convection_x << "_cy_"
                 << convection_y << "_ay_" << diffusion_y << "_" << argv[6];
    stringstream coefffilename;
    coefffilename << "s_adwav_coeff_realline_pde2d_" << argv[1] << "_" << argv[2] << "_"
                  << argv[3] << "_" << argv[4] << "_" << argv[5] << "_cx_" << convection_x << "_cy_"
                  << convection_y << "_ay_" << diffusion_y << "_" << argv[6] << ".dat";

    //Righthand side construction for tensor solution
    if (strcmp(argv[1],"CDF")==0) {
        cout << "Not implemented yet" << endl;
        return 0;
    }
    else if (strcmp(argv[1],"MW")==0) {
        cout << "Not implemented yet" << endl;
        return 0;
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D       SparseMW_basis_x(d,j0_x);
        SparseMW_Basis1D       SparseMW_basis_y(d,j0_y);
        SparseMW_Basis2D       SparseMW_basis2d(SparseMW_basis_x,SparseMW_basis_y);
        SparseMW_MA            SparseMW_A(SparseMW_basis2d,reaction,convection_x,convection_y,diffusion_y);
        //SparseMW_MA            SparseMW_A(SparseMW_basis2d,1.);
        SparseMW_H1_MA         SparseMW_H1_A(SparseMW_basis2d,1.);
        SparseMW_Prec          SparseMW_P(SparseMW_A);
        SparseMW_H1_Prec       SparseMW_H1_P(SparseMW_H1_A);

        if (example==1 || example==2 || example==3) {
            int order = 20;
            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, reaction, convection_x, convection_y, diffusion_y);

            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);
            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_x(SparseMW_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_y(SparseMW_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
            SparseMW_SumOfSeparableRhsIntegral2D SparseMW_rhsintegral2d(SparseMW_rhsintegral_x,SparseMW_rhsintegral_y);
            SparseMW_SumOfSeparableRhs           SparseMW_F(SparseMW_rhsintegral2d,SparseMW_P);

            SeparableFunction2D<T> PP_SepFunc1(refsol.H1_rhs_x, refsol.sing_pts_x,
                                               refsol.exact_y, refsol.sing_pts_y);
            SeparableFunction2D<T> PP_SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                               refsol.H1_rhs_y, refsol.sing_pts_y);
            SparseMW_SeparableRhsIntegral2D      SparseMW_H1_rhsintegral_x(SparseMW_basis2d, PP_SepFunc1, refsol.H1_deltas_x, no_deltas, order);
            SparseMW_SeparableRhsIntegral2D      SparseMW_H1_rhsintegral_y(SparseMW_basis2d, PP_SepFunc2, no_deltas, refsol.H1_deltas_y, order);
            SparseMW_SumOfSeparableRhsIntegral2D SparseMW_H1_rhsintegral2d(SparseMW_H1_rhsintegral_x,SparseMW_H1_rhsintegral_y);
            SparseMW_H1_SumOfSeparableRhs        SparseMW_H1_F(SparseMW_H1_rhsintegral2d,SparseMW_H1_P);

            SparseMW_S_ADWAV_SOLVER_Optim_SeparableRhs SparseMW_s_adwav_solver_optim
                            (SparseMW_A, SparseMW_F, SparseMW_H1_A, SparseMW_H1_F,
                             contraction, threshTol, cgTol, resTol, NumOfIterations, 1, 1e-2);

            SparseMW_A.clear();
            Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
            if (convection_x==0 && convection_y==0) {
                SparseMW_s_adwav_solver_optim.solve(InitialLambda, u, "cg",
                                                   convfilename.str().c_str(), 0, refsol.H1norm());
            }
            else {
                SparseMW_s_adwav_solver_optim.solve(InitialLambda, u, "gmres",
                                                    convfilename.str().c_str(), 0, refsol.H1norm());
            }

            plotScatterCoeff2D(u, SparseMW_basis_x, SparseMW_basis_y, coefffilename.str().c_str());


            stringstream rhsfilename;
            rhsfilename << "rhs_realline_pde2d_" << argv[1] << "_" << argv[2] << "_" << argv[3]
                        << "_" << argv[4] << "_" << argv[5] << "_cx_" << convection_x << "_cy_"
                        << convection_y << "_ay_" << diffusion_y << "_" << argv[6] << ".dat";

            Coefficients<Lexicographical,T,Index2D> f;
            IndexSet<Index2D> Lambda;
            Lambda = supp(u);

            IndexSet<Index2D> Extension;
            Extension = C(Lambda,1.,SparseMW_basis2d);
            Lambda = Lambda + Extension;
            std::cerr << "  Size of enlarged Lambda: " << Lambda.size() << std::endl;
            Extension = C(Lambda,1.,SparseMW_basis2d);
            Lambda = Lambda + Extension;
            std::cerr << "  Size of enlarged Lambda: " << Lambda.size() << std::endl;
            f = SparseMW_F(Lambda);
            std::cerr << "  Finished computation of f with support size " << Lambda.size() << std::endl;
            Coefficients<AbsoluteValue,T,Index2D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=50; ++k) {
                T eta=pow(2.,(T)-k);
                f = SparseMW_F(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index2D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set2d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();

            cout << "Plot of solution started." << endl;
            plot2D(SparseMW_basis2d, u, SparseMW_P, refsol.exact, -10., 10., -10., 10.,
                   pow2i<T>(-3), plotfilename.str().c_str());
            cout << "Plot of solution finished." << endl;




        }
    }
    else {
        std::cerr << "Not yet implemented for " << argv[1] << std::endl;
        return 0;
    }

    return 0;
}

