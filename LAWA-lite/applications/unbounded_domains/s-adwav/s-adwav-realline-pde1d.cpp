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
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;

//Basis definitions
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;
//typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
//typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptivePDEOperatorOptimized1D<T,Primal,R,CDF>            CDF_MA;
//typedef AdaptivePDEOperatorOptimized1D<T,Orthogonal,R,Multi>      MW_MA;
//typedef AdaptivePDEOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, CDF_MA>       CDF_Prec;
//typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, MW_MA>        MW_Prec;
//typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>  SparseMW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, CDF_Basis1D>                                  CDF_RhsIntegral1D;
//typedef RHSWithPeaks1D<T, MW_Basis1D>                                   MW_RhsIntegral1D;
//typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                             SparseMW_RhsIntegral1D;

typedef RHS1D<T, CDF_RhsIntegral1D, CDF_Prec>                           CDF_Rhs;
//typedef RHS1D<T, MW_RhsIntegral1D, MW_Prec>                             MW_Rhs;
//typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                 SparseMW_Rhs;

//Algorithm definition
typedef S_ADWAV<T,Index1D, CDF_Basis1D, CDF_MA, CDF_Rhs>                CDF_S_ADWAV_SOLVER;
//typedef S_ADWAV<T,Index1D, MW_Basis1D, MW_MA, MW_Rhs>                   MW_S_ADWAV_SOLVER;
//typedef S_ADWAV<T,Index1D, SparseMW_Basis1D, SparseMW_MA, SparseMW_Rhs> SparseMW_S_ADWAV_SOLVER;

int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0; bool w_XBSpline;
    int example=atoi(argv[5]);
    int NumOfIterations=atoi(argv[6]);

    T reaction = 1., convection = 15., diffusion = 1.;
    T contraction = 1.;
    T threshTol = 0.1, cgTol = 0.1*threshTol, resTol=1e-4;

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,reaction,0.,diffusion);
    Function<T>        rhsFct(refsol.rhs,refsol.sing_pts);

    if      (strcmp(argv[4],"-inf")==0) {
        std::cerr << "No minimal level implemented for S-ADWAV PDE problem. Exit." << std::endl;
        exit(1);

    }
    else {
        j0=atoi(argv[4]);
    }

    if ((JMINOFFSET + j0) < 0) {
        cout << "Please re-adjust #JMINOFFSET in index.h" << endl;
        exit(1);
    }

    stringstream convfilename;
    convfilename << "s_adwav_conv_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << j0 <<  "_" << reaction << "_" << convection << "_"
                 << diffusion << "_" << argv[5] << ".dat";

    stringstream plotfilename;
    plotfilename << "s_adwav_plot_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << j0 <<  "_" << reaction  << "_" << convection << "_"
                 << diffusion << "_" << argv[5] << ".dat";


    cout << "Initializing S-ADWAV, jmin = " << j0 << endl;

    IndexSet<Index1D> InitialLambda;
    InitialLambda.insert(Index1D(j0,1,XBSpline));

    CDF_Basis1D             CDF_basis(d,d_,j0);
    CDF_MA                  CDF_A(CDF_basis, reaction, convection, diffusion);
    CDF_Prec                CDF_prec(CDF_A);
    CDF_RhsIntegral1D       CDF_rhsintegral1d(CDF_basis, rhsFct, refsol.deltas, 40);
    CDF_Rhs                 CDF_F(CDF_rhsintegral1d,CDF_prec);
    CDF_S_ADWAV_SOLVER      CDF_s_adwav_solver(CDF_basis, CDF_A, CDF_F,
                                               contraction, threshTol, cgTol, resTol,
                                               NumOfIterations, 2, 1e-7, 100000);

    CDF_s_adwav_solver.solve(InitialLambda, "gmres", convfilename.str().c_str(),
                             2, refsol.H1norm());

    cout << "Plot of solution started." << endl;
    plot<T, CDF_Basis1D, CDF_Prec>(CDF_basis, CDF_s_adwav_solver.solutions[NumOfIterations-1], CDF_prec, refsol.u, refsol.d_u, -30., 120.,
                                   0.01, plotfilename.str().c_str());
    cout << "Plot of solution finished." << endl;

    return 0;


}
