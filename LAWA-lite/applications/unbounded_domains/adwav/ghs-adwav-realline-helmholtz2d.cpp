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
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>

#define ROW_SIZE 4*8192
#define COL_SIZE 4*2048

typedef long double T;
using namespace lawa;
using namespace std;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

///  Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                               const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator         const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator           const_coeff2d_abs_it;

///  Basis definitions
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef TensorBasis2D<Adaptive, MW_Basis1D, MW_Basis1D>                 MW_Basis2D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;
typedef TensorBasis2D<Adaptive, SparseMW_Basis1D,SparseMW_Basis1D>      SparseMW_Basis2D;

///  Operator definition
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,R,Multi,
                                               Orthogonal,R,Multi>      MW_MA;
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,R,SparseMulti,
                                               Primal,R,SparseMulti>    SparseMW_MA;

///  (Adaptive) preconditioner definition (allowing for storing computed values)
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D, MW_MA>        MW_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_MA>   SparseMW_Prec;

///  Righthandsides definitions for separable right-hand sides ($f(x_1,x_2) = f(x_1) \cdot f(x_2)$)
///  or sums of such functions
typedef SeparableRHS2D<T,MW_Basis2D>                                    MW_SeparableRhsIntegral2D;
typedef SeparableRHS2D<T,SparseMW_Basis2D>                              SparseMW_SeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,MW_SeparableRhsIntegral2D,
                             MW_SeparableRhsIntegral2D>                 MW_SumOfSeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,SparseMW_SeparableRhsIntegral2D,
                             SparseMW_SeparableRhsIntegral2D>           SparseMW_SumOfSeparableRhsIntegral2D;

typedef RHS2D<T,MW_SumOfSeparableRhsIntegral2D,MW_Prec>                 MW_SumOfSeparableRhs;
typedef RHS2D<T,SparseMW_SumOfSeparableRhsIntegral2D,SparseMW_Prec>     SparseMW_SumOfSeparableRhs;

///  Righthandsides definitions for non-separable right-hand side functions with different quadrature
///  ruls (full and sparse tensor product rules).
typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, SparseGridGP>         MW_NonSeparableRhsIntegralSG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, SparseGridGP>   SparseMW_NonSeparableRhsIntegralSG2D;

typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, FullGridGL>           MW_NonSeparableRhsIntegralFG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, FullGridGL>     SparseMW_NonSeparableRhsIntegralFG2D;

typedef SumOfThreeRHSIntegrals<T, Index2D,
                               MW_NonSeparableRhsIntegralFG2D>          MW_SumOfNonSeparableRhsIntegral2D;
typedef SumOfThreeRHSIntegrals<T, Index2D,
                               SparseMW_NonSeparableRhsIntegralFG2D>    SparseMW_SumOfNonSeparableRhsIntegral2D;

typedef RHS2D<T,MW_NonSeparableRhsIntegralSG2D, MW_Prec>                MW_NonSeparableRhs2D;
typedef RHS2D<T,SparseMW_NonSeparableRhsIntegralSG2D, SparseMW_Prec>    SparseMW_NonSeparableRhs2D;

typedef RHS2D<T,MW_SumOfNonSeparableRhsIntegral2D,MW_Prec>              MW_SumOfNonSeparableRhs2D;
typedef RHS2D<T,SparseMW_SumOfNonSeparableRhsIntegral2D,SparseMW_Prec>  SparseMW_SumOfNonSeparableRhs2D;


///  Definition of APPLY-AWGM solvers depending on the wavelet basis discretization and the right-hand
///  side type
typedef GHS_ADWAV<T, Index2D, MW_MA, MW_SumOfSeparableRhs>              MW_GHS_ADWAV_SOLVER_SumofSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA, SparseMW_SumOfSeparableRhs>  SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs;
typedef GHS_ADWAV<T, Index2D, MW_MA, MW_NonSeparableRhs2D>              MW_GHS_ADWAV_SOLVER_NonSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA, SparseMW_NonSeparableRhs2D>  SparseMW_GHS_ADWAV_SOLVER_NonSeparableRhs;
typedef GHS_ADWAV<T, Index2D, MW_MA, MW_SumOfNonSeparableRhs2D>         MW_GHS_ADWAV_SOLVER_SumNonSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA,
                  SparseMW_SumOfNonSeparableRhs2D>                      SparseMW_GHS_ADWAV_SOLVER_SumNonSeparableRhs;


int main (int argc, char *argv[]) {
    if (argc!=8) {
        cout << "usage " << argv[0] << " basistype d d_ j0_x j0_y example NumOfIterations" << endl; exit(1);
    }
    cout.precision(20);

    ///  Wavelet basis parameters
    int d   =atoi(argv[2]);
    int d_  =atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);

    if (((JMINOFFSET + j0_x) < 0) || ((JMINOFFSET + j0_y) < 0)) {
        cout << "Please re-adjust #JMINOFFSET in index.h" << endl;
        exit(1);
    }


    T c = 1.;
    T diffusion_y = 1.;

    ///  Number of the reference example to be considered
    int example=atoi(argv[6]);

    ///  Maximum number of iterations of the APPLY-AWGM.
    int NumOfIterations=atoi(argv[7]);

    ///  File name for the required rhs file
    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_pde2d_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_cx_0_cy_0_ay_" << diffusion_y << "_"
                << argv[6] << ".dat";

    /*
    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";
    */
    stringstream convfilename;
    convfilename << "ghs_adwav_conv_realline_helmholtz2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_c_" << c << "_ay_"
                 << diffusion_y << "_" << argv[6] << ".dat";

    ///  Integration order (to be adjusted in the sequel)
    int order=20;

    if (strcmp(argv[1],"MW")==0) {

        ///  Initialization of wavelet basis and operators
        MW_Basis1D MW_basis_x(d,j0_x);
        MW_Basis1D MW_basis_y(d,j0_y);
        MW_Basis2D MW_basis2d(MW_basis_x,MW_basis_y);
        MW_MA      MW_A(MW_basis2d,1.);
        MW_Prec    MW_P(MW_A);

        int assemble_matrix = 2;    //assemble matrix, use operator internal method

        if (example==1 || example==2 || example==3) {

            ///  Setting up the right-hand side for separable reference solution
            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1.,0., 0., 1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
            MW_SeparableRhsIntegral2D MW_rhsintegral_x(MW_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
            MW_SeparableRhsIntegral2D MW_rhsintegral_y(MW_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
            MW_SumOfSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_x,MW_rhsintegral_y);
            MW_SumOfSeparableRhs MW_F(MW_rhsintegral2d,MW_P);
            MW_F.readIndexSets(rhsfilename.str().c_str());

            ///  Initialization of the APPLY-AWGM solver
            MW_GHS_ADWAV_SOLVER_SumofSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);

            ///  Calling the APPLY-AWGM solver
            Coefficients<Lexicographical,T,Index2D> u;
            u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                          NumOfIterations, refsol.H1norm());

        }
        else if (example==4) {

            ///  Setting up the right-hand side for separable reference solution
            int order = 31;
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(1, 1.);
            Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
            MW_NonSeparableRhsIntegralSG2D MW_rhsintegral2d(MW_basis2d, Func2d, order);
            MW_NonSeparableRhs2D           MW_F(MW_rhsintegral2d,MW_P);
            MW_F.readIndexSets(rhsfilename.str().c_str());

            ///  Initialization of the APPLY-AWGM solver
            MW_GHS_ADWAV_SOLVER_NonSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);

            ///  Calling the APPLY-AWGM solver
            Coefficients<Lexicographical,T,Index2D> u;
            u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                          NumOfIterations, refsol.H1norm());

        }
        else if (example==5) {
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(2, 1.);
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);

            if (d==2 || d==3) {
                int order = 40; // required for convergence on high levels -> Attention: slow!
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_reaction(MW_basis2d, Func2d, order);
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_diffusion_x(MW_basis2d, Func2d_x, order, 1, 0);
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_diffusion_y(MW_basis2d, Func2d_y, order, 0, 1);
                MW_SumOfNonSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_diffusion_x,
                                                                   MW_rhsintegral_diffusion_y,
                                                                   MW_rhsintegral_reaction);
                MW_SumOfNonSeparableRhs2D         MW_F(MW_rhsintegral2d,MW_P);
                MW_F.readIndexSets(rhsfilename.str().c_str());

                MW_GHS_ADWAV_SOLVER_SumNonSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);

                Coefficients<Lexicographical,T,Index2D> u;
                u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                              NumOfIterations, refsol.H1norm());
            }
        }
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D SparseMW_basis_x(d,j0_x);
        SparseMW_Basis1D SparseMW_basis_y(d,j0_y);
        SparseMW_Basis2D SparseMW_basis2d(SparseMW_basis_x,SparseMW_basis_y);
        SparseMW_MA      SparseMW_A(SparseMW_basis2d,1.,diffusion_y);
        SparseMW_Prec    SparseMW_P(SparseMW_A);
        int assemble_matrix = 0;    //compute matrix vector product solely on index set
        if (example==1 || example==2 || example==3) {

            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1., 0., 0., diffusion_y);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_x(SparseMW_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_y(SparseMW_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
            SparseMW_SumOfSeparableRhsIntegral2D SparseMW_rhsintegral2d(SparseMW_rhsintegral_x,SparseMW_rhsintegral_y);
            SparseMW_SumOfSeparableRhs           SparseMW_F(SparseMW_rhsintegral2d,SparseMW_P);
            SparseMW_F.readIndexSets(rhsfilename.str().c_str());

            SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs SparseMW_ghs_adwav_solver(SparseMW_A,SparseMW_F,true,assemble_matrix);
            T alpha = 0.6, omega = 0.2, gamma = 0.15, theta = 2*omega/(1+omega);
            SparseMW_ghs_adwav_solver.setParameters(alpha, omega, gamma, theta);

            Coefficients<Lexicographical,T,Index2D> u;
            u = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                                NumOfIterations, refsol.Energynorm());
        }
    }
    return 0;
}

