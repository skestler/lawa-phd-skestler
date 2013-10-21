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

typedef double T;
using namespace lawa;
using namespace std;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                               const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator         const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator           const_coeff2d_abs_it;

//Basis definitions
typedef Basis<T,Primal,Interval,SparseMulti>                            SparseMW_IntervalBasis1D;
typedef Basis<T,Primal,RPlus,SparseMulti>                               SparseMW_RPlusBasis1D;
typedef TensorBasis2D<Adaptive, SparseMW_IntervalBasis1D,
                                SparseMW_RPlusBasis1D>                  SparseMW_Basis2D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Interval,SparseMulti,
                                        Primal,RPlus,SparseMulti>       SparseMW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_MA>   SparseMW_Prec;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,SparseMW_Basis2D>                              SparseMW_SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SparseMW_SeparableRhsIntegral2D,
                             SparseMW_SeparableRhsIntegral2D>           SparseMW_SumOfSeparableRhsIntegral2D;

typedef RHS2D<T,SparseMW_SumOfSeparableRhsIntegral2D,SparseMW_Prec>     SparseMW_SumOfSeparableRhs;

//Algorithm definition
typedef GHS_ADWAV<T, Index2D, SparseMW_MA, SparseMW_SumOfSeparableRhs>  SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs;


int main (int argc, char *argv[]) {
    if (argc!=8) {
        cout << "usage " << argv[0] << " basistype d d_ j0_x j0_y example NumOfIterations" << endl; exit(1);
    }
    cout.precision(20);

    int d   =atoi(argv[2]);
    int d_  =atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);
    T reaction=1., convection_x=0., convection_y=0., diffusion_y=1.;
    int example=atoi(argv[6]);
    int NumOfIterations=atoi(argv[7]);

    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_interval_rplus_pde2d_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_cx_" << convection_x << "_cy_"
                << convection_y << "_ay_" << diffusion_y << "_" << argv[6] << ".dat";

    stringstream convfilename;
    convfilename << "ghs_adwav_conv_interval_rplus_helmholtz2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_c_" << reaction << "_ay_"
                 << diffusion_y << "_" << argv[6] << ".dat";

    stringstream coefffilename;
    coefffilename << "ghs_adwav_coeff_interval_rplus_helmholtz2d_" << argv[1] << "_" << argv[2] << "_"
                  << argv[3] << "_" << argv[4] << "_" << argv[5] << "_c_" << reaction << "_ay_"
                  << diffusion_y << "_" << argv[6] << ".dat";

    int order=40;

    SparseMW_IntervalBasis1D SparseMW_basis_x(d,j0_x);
    SparseMW_basis_x.enforceBoundaryCondition<DirichletBC>();
    SparseMW_RPlusBasis1D    SparseMW_basis_y(d,j0_y);
    SparseMW_basis_y.enforceBoundaryCondition<DirichletBC>();
    SparseMW_Basis2D SparseMW_basis2d(SparseMW_basis_x,SparseMW_basis_y);
    SparseMW_MA      SparseMW_A(SparseMW_basis2d,reaction,diffusion_y);
    SparseMW_Prec    SparseMW_P(SparseMW_A);

    TensorRefSols_PDE_Interval_RPlus<T> refsol;
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
    SparseMW_F.readIndexSets(rhsfilename.str().c_str());

    int assemble_matrix = 0;    //compute matrix vector soley on indexset
    SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs SparseMW_ghs_adwav_solver(SparseMW_A,SparseMW_F,true,assemble_matrix);
    T alpha = 0.6, omega = 0.2, gamma = 0.15, theta = 2*omega/(1+omega);
    SparseMW_ghs_adwav_solver.setParameters(alpha, omega, gamma, theta);

    Coefficients<Lexicographical,T,Index2D> u;
    u = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                        NumOfIterations, refsol.Energynorm());

    plot2D<T,SparseMW_Basis2D,SparseMW_Prec>(SparseMW_basis2d, u, SparseMW_P, refsol.exact,
                               0., 1., 0., 10., pow2i<T>(-5), "example_interval_rplus");

    return 0;
}

