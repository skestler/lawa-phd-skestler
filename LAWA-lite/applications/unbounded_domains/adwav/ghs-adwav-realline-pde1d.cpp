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
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptivePDEOperatorOptimized1D<T,Primal,R,SparseMulti>          SparseMW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_PP_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>  SparseMW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                             SparseMW_RhsIntegral1D;

typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                 SparseMW_Rhs;

typedef GHS_NONSYM_ADWAV<T, Index1D, SparseMW_MA, SparseMW_Rhs,
                         SparseMW_PP_MA, SparseMW_Rhs>               SparseMW_GHS_NONSYM_ADWAV_SOLVER;


int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0=atoi(argv[4]);

    T reaction = 1., convection = 10., diffusion = 1.;
    int example=atoi(argv[5]);
    int NumOfIterations=atoi(argv[6]);

    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";
    cout << rhsfilename.str().c_str() << endl;
    stringstream convfilename;
    convfilename << "ghs_adwav_conv_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";
    stringstream plotfilename;
    plotfilename << "ghs_adwav_plot_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,reaction,convection,diffusion);
    Function<T>                 rhsFct(refsol.rhs,refsol.sing_pts);
    Function<T>                 PP_rhsFct(refsol.H1_rhs,refsol.sing_pts);

    SparseMW_Basis1D            SparseMW_basis(d,j0);
    SparseMW_MA                 SparseMW_A(SparseMW_basis,reaction,convection,diffusion);
    SparseMW_PP_MA              SparseMW_PP_A(SparseMW_basis,1.);
    SparseMW_Prec               SparseMW_prec(SparseMW_A);
    SparseMW_RhsIntegral1D      SparseMW_rhsintegral1d(SparseMW_basis, rhsFct, refsol.deltas, 17);
    SparseMW_RhsIntegral1D      SparseMW_PP_rhsintegral1d(SparseMW_basis, PP_rhsFct, refsol.H1_deltas, 17);
    SparseMW_Rhs                SparseMW_F(SparseMW_rhsintegral1d,SparseMW_prec);
    SparseMW_Rhs                SparseMW_PP_F(SparseMW_PP_rhsintegral1d,SparseMW_prec);
    SparseMW_GHS_NONSYM_ADWAV_SOLVER   SparseMW_ghs_adwav_solver(SparseMW_A, SparseMW_F,
                                                                 SparseMW_PP_A, SparseMW_PP_F,true);


    if (SparseMW_F.readIndexSets(rhsfilename.str().c_str()) ) {
        cout << "Index sets for rhs read... Ready to start."  << endl;
    }
    else {
        cout << "RHS: Could not open file." << endl;
        return 0;
    }

    Coefficients<Lexicographical,T,Index1D> w;
    w = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-10, convfilename.str().c_str(),
                                        NumOfIterations, refsol.H1norm());
    cout << "Plot of solution started." << endl;
    plot<T, SparseMW_Basis1D, SparseMW_Prec>(SparseMW_basis, w, SparseMW_prec, refsol.u,
                                             refsol.d_u, -80., 80., pow2i<T>(-5),
                                             plotfilename.str().c_str());
    cout << "Plot of solution finished." << endl;

    return 0;
}
