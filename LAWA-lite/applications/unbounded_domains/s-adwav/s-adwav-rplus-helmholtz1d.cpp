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

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >            SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >          DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                                DenseVectorT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                                   const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator             const_coeff1d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator               const_coeff1d_abs_it;

//Basis definitions
typedef Basis<T,Primal,RPlus,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,RPlus,SparseMulti>    SparseMW_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>      SparseMW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                                 SparseMW_RhsIntegral1D;

typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                     SparseMW_Rhs;

//Algorithm definition
typedef S_ADWAV<T,Index1D, SparseMW_Basis1D, SparseMW_MA, SparseMW_Rhs>     SparseMW_S_ADWAV_SOLVER;

int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0;
    int example=atoi(argv[5]);
    int NumOfIterations=atoi(argv[6]);

    T c = 1.;
    T contraction = 1.;
    T threshTol = 0.1, cgTol = 0.1*threshTol, resTol=1e-4;

    RefSols_PDE_RPlus1D<T> refsol;
    refsol.setExample(example,1.,0.,c);
    Function<T>        rhsFct(refsol.rhs,refsol.sing_pts);

    j0=atoi(argv[4]);

    if ((JMINOFFSET + j0) < 0) {
        cout << "Please re-adjust #JMINOFFSET in index.h" << endl;
        exit(1);
    }

    stringstream convfilename;
    convfilename << "s_adwav_conv_rplus_helmholtz1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << j0 <<  "_" << c << "_" << argv[5] << ".dat";

    stringstream plotfilename;
    plotfilename << "s_adwav_plot_rplus_helmholtz1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << j0 <<  "_" << c << "_" << argv[5] << ".dat";


    cout << "Initializing S-ADWAV, jmin = " << j0 << endl;

    IndexSet<Index1D> InitialLambda;
    InitialLambda.insert(Index1D(j0,1,XBSpline));

    if (strcmp(argv[1],"CDF")==0) {
        cerr << "Not implemented yet." << endl;
        exit(1);
    }
    else if (strcmp(argv[1],"MW")==0) {
        cerr << "Not implemented yet." << endl;
        exit(1);
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D            SparseMW_basis(d,j0);
        SparseMW_basis.enforceBoundaryCondition<DirichletBC>();
        SparseMW_MA                 SparseMW_A(SparseMW_basis,c);
        SparseMW_Prec               SparseMW_prec(SparseMW_A);
        SparseMW_RhsIntegral1D      SparseMW_rhsintegral1d(SparseMW_basis, rhsFct, refsol.deltas, 40);
        SparseMW_Rhs                SparseMW_F(SparseMW_rhsintegral1d,SparseMW_prec);
        SparseMW_S_ADWAV_SOLVER     SparseMW_s_adwav_solver(SparseMW_basis, SparseMW_A, SparseMW_F,
                                                            contraction, threshTol, cgTol, resTol,
                                                            NumOfIterations, 2, 1e-7, 100000);


        SparseMW_s_adwav_solver.solve(InitialLambda, "cg",convfilename.str().c_str(),
                                      2, refsol.H1norm());

        cout << "Plot of solution started." << endl;
        plot<T, SparseMW_Basis1D, SparseMW_Prec>(SparseMW_basis, SparseMW_s_adwav_solver.solutions[NumOfIterations-1], SparseMW_prec, refsol.u, refsol.d_u, 0., 80.,
                                       0.01, plotfilename.str().c_str());
        cout << "Plot of solution finished." << endl;
    }

    return 0;

}
