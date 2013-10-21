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

///  Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;

///  Basis definitions
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

///  (Adaptive) operator definitions (allowing for storing computed values)
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>            CDF_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,R,Multi>      MW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_MA;

///  (Adaptive) preconditioner definitions (allowing for storing computed values)
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, CDF_MA>       CDF_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, MW_MA>        MW_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>  SparseMW_Prec;

///  Right-hand side integral definitions
typedef RHSWithPeaks1D<T, CDF_Basis1D>                                  CDF_RhsIntegral1D;
typedef RHSWithPeaks1D_WO_XBSpline<T>                                   CDF_RhsIntegral1D_WO_XBSpline;
typedef RHSWithPeaks1D<T, MW_Basis1D>                                   MW_RhsIntegral1D;
typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                             SparseMW_RhsIntegral1D;

///  Right-hand side definitions
typedef RHS1D<T, CDF_RhsIntegral1D, CDF_Prec>                           CDF_Rhs;
typedef RHS1D<T, CDF_RhsIntegral1D_WO_XBSpline, CDF_Prec>               CDF_Rhs_WO_XBSpline;
typedef RHS1D<T, MW_RhsIntegral1D, MW_Prec>                             MW_Rhs;
typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                 SparseMW_Rhs;

///  Definition of solvers (as for all other definitions before, depending on the wavelet basis).
///  In particular, we need to distinguish between the cases where the minimal level $j_0 \neq -\infty$
///  and where we do not have a minimal level.
typedef GHS_ADWAV<T, Index1D, CDF_MA, CDF_Rhs>                          CDF_GHS_ADWAV_SOLVER;
typedef GHS_ADWAV<T, Index1D, CDF_MA, CDF_Rhs_WO_XBSpline>              CDF_GHS_ADWAV_SOLVER_WO_XBSpline;
typedef GHS_ADWAV<T, Index1D, MW_MA, MW_Rhs>                            MW_GHS_ADWAV_SOLVER;
typedef GHS_ADWAV<T, Index1D, SparseMW_MA, SparseMW_Rhs>                SparseMW_GHS_ADWAV_SOLVER;

int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example max_its" << endl; exit(1);
    }
    cout.precision(3);

    ///  Wavelet parameters
    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0; bool w_XBSpline;

    ///  Constant appearing the operator $-\Delta + c\cdot \textrm{Id}$.
    T c = 1.;

    ///  Number of reference example and maximum number of AWGM iterations
    int example=atoi(argv[5]);
    int NumOfIterations=atoi(argv[6]);

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,c);
    Function<T>        rhsFct(refsol.rhs,refsol.sing_pts);

    ///  Definition of required rhs file name required for $\textbf{RHS}$. Moreover, we distunguish
    ///  the case between minimal level $j_0 \neq -\infty$ and no minimal level
    stringstream rhsfilename;
//    chdir("./rhs");
    if      (strcmp(argv[4],"-inf")==0) {
        j0=0;             w_XBSpline=false;
        rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << argv[4] << "_" << c << "_" << argv[5] << ".dat";
    }
    else if (strcmp(argv[4],"best")==0) {
        j0=refsol.getMinimalLevel(d,d_); w_XBSpline=true;
        rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << j0 << "_" << c << "_" << argv[5] << ".dat";
    }
    else {
        j0=atoi(argv[4]); w_XBSpline=true;
        rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << argv[4] << "_" << c << "_" << argv[5] << ".dat";
    }

    stringstream plotfilename;
    plotfilename << "ghs_adwav_plot_realline_helmholtz1d_" << argv[1] << "_" << argv[2] << "_"
                  << argv[3] << "_" << argv[4] << "_" << c << "_"
                  << argv[5] << "_" << argv[6] << ".dat";

    stringstream convfilename;
    convfilename << "ghs_adwav_conv_realline_helmholtz1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << c << "_" << argv[5] << ".dat";

    Coefficients<Lexicographical,T,Index1D> u;

    if (strcmp(argv[1],"CDF")==0) {

        ///  Initialization of basis and operators
        CDF_Basis1D             CDF_basis(d,d_,j0);
        CDF_MA                  CDF_A(CDF_basis,w_XBSpline,c);
        CDF_Prec                CDF_prec(CDF_A);

        bool optimizedGrow =   true;
        int assembleMatrix  =   2;

        if (w_XBSpline) {
            ///  Initialization of rhs when $j_0 \neq -\infty$
            CDF_RhsIntegral1D       CDF_rhsintegral1d(CDF_basis, rhsFct, refsol.deltas, 70);
            CDF_Rhs                 CDF_F(CDF_rhsintegral1d,CDF_prec);

            ///  Initialization of the AWGM solver when $j_0 \neq -\infty$
            CDF_GHS_ADWAV_SOLVER    CDF_ghs_adwav_solver(CDF_A, CDF_F,optimizedGrow,assembleMatrix);

            if (CDF_F.readIndexSets(rhsfilename.str().c_str()) ) {
                cout << "Index sets for rhs read... Ready to start."  << endl;
            }
            else {
                cout << "RHS: Could not open file." << endl;
                return 0;
            }

            ///  Calling the AWGM solver when $j_0 \neq -\infty$
            u = CDF_ghs_adwav_solver.SOLVE(CDF_F.norm_estimate, 1e-10,  convfilename.str().c_str(),
                                           NumOfIterations, refsol.H1norm());

        }
        else {
            T left_bound = 0., right_bound = 0.;
            ///  Initialization of rhs when $j_0 = -\infty$
            int J_plus_smooth = 0, J_minus_smooth = 0, J_plus_singular = 0, J_minus_singular = 0;
            bool singular_integral=false;
            refsol.getRHS_WO_XBSplineParameters(d, d_, left_bound, right_bound,
                                                J_plus_smooth, J_minus_smooth,
                                                J_plus_singular, J_minus_singular, singular_integral);
            CDF_RhsIntegral1D_WO_XBSpline       CDF_rhsintegral1d_WO_XBSpline
                                                  (CDF_basis.psi, refsol.rhs, refsol.sing_pts,
                                                   refsol.deltas, left_bound, right_bound,
                                                   4., 20);
            CDF_Rhs_WO_XBSpline                 CDF_F_WO_XBSpline(CDF_rhsintegral1d_WO_XBSpline,
                                                                  CDF_prec);

            ///  Initialization of the AWGM solver when $j_0 = -\infty$
            CDF_GHS_ADWAV_SOLVER_WO_XBSpline    CDF_ghs_adwav_solver_WO_XBSpline(CDF_A,CDF_F_WO_XBSpline,optimizedGrow,
                                                                                 assembleMatrix);

            if (CDF_F_WO_XBSpline.readIndexSets(rhsfilename.str().c_str()) ) {
                cout << "Index sets for rhs read... Ready to start."  << endl;
            }
            else {
                cout << "RHS: Could not open file." << endl;
                return 0;
            }

            ///  Calling the AWGM solver when $j_0 \neq -\infty$
            u = CDF_ghs_adwav_solver_WO_XBSpline.SOLVE(CDF_F_WO_XBSpline.norm_estimate, 1e-10,
                                                       convfilename.str().c_str(),
                                                       NumOfIterations, refsol.H1norm());

        }

        cout << "Plot of solution started." << endl;
        plot<T, CDF_Basis1D, CDF_Prec>(CDF_basis, u, CDF_prec, refsol.u, refsol.d_u, -80., 80.,
                                       pow2i<T>(-5), plotfilename.str().c_str());
        cout << "Plot of solution finished." << endl;
    }
    else if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D              MW_basis(d,j0);
        MW_MA                   MW_A(MW_basis,c);
        MW_Prec                 MW_prec(MW_A);
        MW_RhsIntegral1D        MW_rhsintegral1d(MW_basis, rhsFct, refsol.deltas, 20);
        MW_Rhs                  MW_F(MW_rhsintegral1d,MW_prec);
        bool optimizedGrow =    true;
        int assembleMatrix  =   2;
        MW_GHS_ADWAV_SOLVER     MW_ghs_adwav_solver(MW_A, MW_F, optimizedGrow, assembleMatrix);

        if (MW_F.readIndexSets(rhsfilename.str().c_str()) ) {
            cout << "Index sets for rhs read... Ready to start."  << endl;
        }
        else {
            cout << "RHS: Could not open file." << endl;
            return 0;
        }
        u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-10,  convfilename.str().c_str(),
                                      NumOfIterations, refsol.H1norm());
        cout << "Plot of solution started." << endl;
        plot<T, MW_Basis1D, MW_Prec>(MW_basis, u, MW_prec, refsol.u,
                                     refsol.d_u, -20., 20., pow2i<T>(-5),
                                     plotfilename.str().c_str());
        cout << "Plot of solution finished." << endl;

    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D            SparseMW_basis(d,j0);
        SparseMW_MA                 SparseMW_A(SparseMW_basis,c);
        SparseMW_Prec               SparseMW_prec(SparseMW_A);
        SparseMW_RhsIntegral1D      SparseMW_rhsintegral1d(SparseMW_basis, rhsFct, refsol.deltas, 10);
        SparseMW_Rhs                SparseMW_F(SparseMW_rhsintegral1d,SparseMW_prec);
        bool optimizedGrow =    true;
        int assembleMatrix  =   2;
        SparseMW_GHS_ADWAV_SOLVER   SparseMW_ghs_adwav_solver(SparseMW_A, SparseMW_F, optimizedGrow, assembleMatrix);

        T alpha = 0.6, omega = 0.2, gamma = 0.15, theta = 2*omega/(1+omega);
        SparseMW_ghs_adwav_solver.setParameters(alpha, omega, gamma, theta);
        if (SparseMW_F.readIndexSets(rhsfilename.str().c_str()) ) {
            cout << "Index sets for rhs read... Ready to start."  << endl;
        }
        else {
            cout << "RHS: Could not open file." << endl;
            return 0;
        }
        u = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                            NumOfIterations, refsol.H1norm());
        cout << "Plot of solution started." << endl;
        plot<T, SparseMW_Basis1D, SparseMW_Prec>(SparseMW_basis, u, SparseMW_prec, refsol.u,
                                                 refsol.d_u, -20., 20., pow2i<T>(-5),
                                                 plotfilename.str().c_str());
        cout << "Plot of solution finished." << endl;
    }

    return 0;
}


