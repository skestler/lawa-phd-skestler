#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

// Several typedefs for notational convenience.

// Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<long double> >              DenseVectorLD;

/// Wavelet basis over an interval: Dijkema interval basis and $L_2$-orthonormal multiwavelet basis
/// are possible.
typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
//typedef Basis<T, Primal, Interval, Dijkema>                       PrimalBasis;

/// Refinement basis containing $\bar \Phi_j$ - the collection of refinement B-Splines
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

typedef Integral<Gauss,PrimalBasis,PrimalBasis>                     MultiWaveletIntegral;
typedef Integral<Gauss,RefinementBasis,RefinementBasis>             MultiRefinementIntegral;

/// Check if refinement coefficients for refinement B-Splines, i.e., $\bar \phi_{\bar j,k} =
/// \sum_{k} \bar a_{k,m}^{\bar j} \bar \phi_{\bar j+1,k}$, are correct.
void
test_refinementOfBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv);

/// Check if refinement coefficients for scaling functions, i.e., $\phi_{j,k} =
/// \sum_{k} \bar c_{k,m}^{j} \bar \phi_{\bar j+j,k}$, are correct.
void
test_refinementOfScaling(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv);

/// Check if refinement coefficients for wavelets, i.e., $\psi_{j,k} =
/// \sum_{k} \bar b_{k,m}^{j} \bar \phi_{\bar j+j+1,k}$, are correct.
/// Observe that we use a higher level $\bar j+j+1$ here as wavelets on level $j$ are (finite)
/// linear combinations of scaling functions on level $j+1$.
void
test_refinementOfWavelet(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv);

/// Check if refinement B-Spline neighbors $\bar \phi_{\bar j,m}$ for a given refinement
/// B-Spline $\bar \phi_{\bar j,k}$ with $|\mathrm{supp}\,\bar \phi_{\bar j,m} \cap
/// \mathrm{supp}\,\bar \phi_{\bar j,k}|>0$ are computed correctly.
void
test_getBSplineNeighborsForBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis);

/// Check if wavelet neighbors $\psi_{j,m}$ for a given refinement B-Spline $\bar \phi_{j+\bar j+1,k}$
/// with $|\mathrm{supp}\,\psi_{j,m} \cap \mathrm{supp}\,\bar \phi_{j+\bar j+1,k}|>0$ are
/// computed correctly. Here, $j+\bar j+1$ should be the smallest level satisfying
/// $\mathrm{span}\, \bar \Phi_{j+\bar j+1} \supset \mathrm{span}\, \Psi_j$.
void
test_getWaveletNeighborsForBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis);

/// Check if scaling function neighbors $\phi_{j,m}$ for a given scaling function $\phi_{j,k}$
/// with $|\mathrm{supp}\,\phi_{j,m} \cap \mathrm{supp}\,\phi_{j,k}|>0$ are
/// computed correctly.
void
test_getScalingNeighborsForScaling(const PrimalBasis &basis);

/// Check if wavelet neighbors $\psi_{j,m}$ for a given scaling function $\phi_{j,k}$
/// with $|\mathrm{supp}\,\psi_{j,m} \cap \mathrm{supp}\,\phi_{j,k}|>0$ are
/// computed correctly.
void
test_getWaveletNeighborsForScaling(const PrimalBasis &basis);

/// Check if refinement B-Spline neighbors $\bar \phi_{j+\bar j+1,m}$ for a given wavelet $\psi_{j,k}$
/// with $|\mathrm{supp}\,\bar \phi_{j+\bar j+1,m} \cap \mathrm{supp}\,\psi_{j,k}|>0$ are
/// computed correctly. Here, $j+\bar j+1$ should be the smallest level satisfying
/// $\mathrm{span}\, \bar \Phi_{j+\bar j+1} \supset \mathrm{span}\, \Psi_j$.
void
test_getBSplineNeighborsForWavelet(const PrimalBasis &basis, const RefinementBasis &refinementbasis);

/// Check if scaling functions neighbors $\phi_{j,m}$ for a given wavelet $\psi_{j,k}$
/// with $|\mathrm{supp}\,\phi_{j,m} \cap \mathrm{supp}\,\psi_{j,k}|>0$ are
/// computed correctly.
void
test_getScalingNeighborsForWavelet(const PrimalBasis &basis);

/// Check if wavelet neighbors $\psi_{j,m}$ for a given wavelet $\psi_{j,k}$
/// with $|\mathrm{supp}\,\psi_{j,m} \cap \mathrm{supp}\,\psi_{j,k}|>0$ are
/// computed correctly.
void
test_getWaveletNeighborsForWavelet(const PrimalBasis &basis);

/// Check if lower wavelet neighbors $\psi_{j-1,m}$ for a given wavelet $\psi_{j,k}$
/// with $|\mathrm{supp}\,\psi_{j-1,m} \cap \mathrm{supp}\,\psi_{j,k}|>0$ are
/// computed correctly.
void
test_getLowerWaveletNeighborsForWavelet(const PrimalBasis &basis);

/// Check if higher wavelet neighbors $\psi_{j+1,m}$ for a given wavelet $\psi_{j,k}$
/// with $|\mathrm{supp}\,\psi_{j+1,m} \cap \mathrm{supp}\,\psi_{j,k}|>0$ are
/// computed correctly.
void
test_getHigherWaveletNeighborsForWavelet(const PrimalBasis &basis);

// Check precision of some generators
void
test_precisionOfGenerators();

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    // wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);
    int J  = atoi(argv[3]);
    int deriv = 0;

    // Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);     // For L2_orthonormal and special MW bases
    //PrimalBasis basis(d, d, j0);     // For biorthogonal wavelet bases
    //if (d>1) basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis& refinementbasis = basis.refinementbasis;

    // Test refinement of B-splines

    //test_refinementOfBSpline(basis, refinementbasis, deriv);

    // Test refinement of multiscaling functions. In case biorthogonal wavelet bases, this is just
    // the same as the test above as here, the scaling function are already B-splines.

    //test_refinementOfScaling(basis, refinementbasis, deriv);

    // Test refinement of multiwavelets: We check the refinement of wavelets in terms of B-splines.

    test_refinementOfWavelet(basis, refinementbasis, deriv);

    //test_getBSplineNeighborsForBSpline(basis, refinementbasis);

    //test_getWaveletNeighborsForBSpline(basis, refinementbasis);

    //test_getScalingNeighborsForScaling(basis);

    //test_getWaveletNeighborsForScaling(basis);


    // The other routines are analog... with different function types instead of B-splines however.

    //test_getBSplineNeighborsForWavelet(basis, refinementbasis);

    //test_getScalingNeighborsForWavelet(basis);

    //test_getWaveletNeighborsForWavelet(basis);

    //test_getLowerWaveletNeighborsForWavelet(basis);

    //test_getHigherWaveletNeighborsForWavelet(basis);


    return 0;
}

void
test_refinementOfBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv)
{
    DenseVectorLD *refCoeffs;
    cout << " ******** Refinement of B-splines **********" << endl;
    for (int j=refinementbasis.mra.j0; j<=refinementbasis.mra.j0+4; ++j) {
        cout << "j = " << j << ": " << refinementbasis.mra.cardI(j) << " " << refinementbasis.mra.rangeI(j) << endl;
        for (int k= refinementbasis.mra.rangeI(j).firstIndex(); k<=refinementbasis.mra.rangeI(j).lastIndex(); ++k) {
            cout << "  k = " << k  << " " << refinementbasis.generator(XBSpline).support(j,k) << " "
                 << refinementbasis.generator(XBSpline).singularSupport(j,k) << endl;
            ofstream plotfile_scaling("refinement_interval.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            for (T x=0.; x<=1.; x+=pow2i<T>(-10-j)) {
                T reference_value = refinementbasis.generator(XBSpline).operator()(x,j,k,deriv);
                T refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                T tmp1 = fabs(reference_value-refinement_value)/fabs(reference_value);
                T tmp2 = fabs(reference_value-refinement_value);
                if (rel_error<=tmp1 && fabs(reference_value)>1e-15) { rel_error = tmp1; x_rel_crit = x;    }
                if (abs_error<=tmp2) { abs_error = tmp2; x_abs_crit = x;    }
                plotfile_scaling << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_scaling.close();
            cout << "Relative error: " << rel_error << ", x_rel_crit = " << x_rel_crit << endl;
            cout << "Absolute error: " << abs_error << ", x_abs_crit = " << x_abs_crit << endl;
            cout << "Please hit enter." << endl;
            getchar();
        }
    }
    cout << " *******************************************" << endl << endl;
}

void
test_refinementOfScaling(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv)
{
    DenseVectorLD *refCoeffs;
    cout << " ******* Refinement of multiscalings *******" << endl;
    for (int j=basis.j0; j<=basis.j0+4; ++j) {
        for (int k=basis.mra.rangeI(j).firstIndex(); k<=basis.mra.rangeI(j).lastIndex(); ++k) {
            ofstream plotfile_scaling("refinement_interval_multiscaling.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = basis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = "
                           << refinement_j << ", refinement_k_first = " << refinement_k_first << endl;
            cout << refinementbasis.mra.rangeI(refinement_j) << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-8-j)) {
                T reference_value = basis.generator(XBSpline)(x,j,k,deriv);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * basis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                if (x==0.5) {
                    cout << "x = 0.5: " << reference_value << " " << refinement_value << endl;
                }
                T tmp1 = fabs(reference_value-refinement_value)/fabs(reference_value);
                T tmp2 = fabs(reference_value-refinement_value);
                if (rel_error<=tmp1 && fabs(reference_value)>1e-15) { rel_error = tmp1; x_rel_crit = x;    }
                if (abs_error<=tmp2) { abs_error = tmp2; x_abs_crit = x;    }
                plotfile_scaling << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_scaling.close();
            cout << "Relative error: " << rel_error << ", x_rel_crit = " << x_rel_crit << endl;
            cout << "Absolute error: " << abs_error << ", x_abs_crit = " << x_abs_crit << endl;
            cout << "Please hit enter." << endl;
            getchar();
        }
    }
    cout << " *******************************************" << endl << endl;
}

void
test_refinementOfWavelet(const PrimalBasis &basis, const RefinementBasis &refinementbasis, int deriv)
{
    DenseVectorLD *refCoeffs;
    cout << " *******************************************" << endl;
    cout << " ******* Refinement of multiwavelets *******" << endl;
    for (int j=basis.j0; j<=basis.j0+10; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            ofstream plotfile_wavelet("refinement_interval_multiwavelet.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = "
                 << refinement_j << ", refinement_k_first = " << refinement_k_first << endl;
            for (T x=0.L; x<=1.L; x+=pow2i<long double>(-8-j)) {
                T reference_value = basis.generator(XWavelet)(x,j,k,deriv);
                T refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {

                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                T tmp1 = fabs(reference_value-refinement_value)/fabs(reference_value);
                T tmp2 = fabs(reference_value-refinement_value);
                if (rel_error<=tmp1 && fabs(reference_value)>1e-15) { rel_error = tmp1; x_rel_crit = x;    }
                if (abs_error<=tmp2) { abs_error = tmp2; x_abs_crit = x;    }
                plotfile_wavelet << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_wavelet.close();
            cout << "Relative error: " << rel_error << ", x_rel_crit = " << x_rel_crit << endl;
            cout << "Absolute error: " << abs_error << ", x_abs_crit = " << x_abs_crit << endl;
            cout << "Please hit enter." << endl;

            if (fabs(abs_error)>1e-12) {
                cout << "Critical error!" << endl;
                getchar();
            }

            //getchar();
        }
    }
    cout << " *******************************************" << endl << endl;
}

void
test_getBSplineNeighborsForBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis)
{
    cout << " ******** BSpline neighbors for BSpline **********" << endl;
    for (int j_bspline1=refinementbasis.j0; j_bspline1<refinementbasis.j0+4; ++j_bspline1) {
        for (long k_bspline1= refinementbasis.mra.rangeI(j_bspline1).firstIndex();
                  k_bspline1<=refinementbasis.mra.rangeI(j_bspline1).lastIndex(); ++k_bspline1) {
            int j_bspline2=0;
            long k_bspline2_first=0L, k_bspline2_last=0L;
            refinementbasis.getBSplineNeighborsForBSpline(j_bspline1, k_bspline1, refinementbasis,
                                                    j_bspline2, k_bspline2_first, k_bspline2_last);

            cout << "BSpline (" << j_bspline1 << "," << k_bspline1 << "): "
                                << j_bspline2 << " , [" << k_bspline2_first << "," << k_bspline2_last << "], "
                                << refinementbasis.mra.rangeI(j_bspline1) << endl;
            for (long k_bspline2=refinementbasis.mra.rangeI(j_bspline2).firstIndex();
                      k_bspline2<k_bspline2_first; ++k_bspline2) {
                if (overlap(refinementbasis.mra.phi.support(j_bspline1,k_bspline1),
                            refinementbasis.mra.phi.support(j_bspline2,k_bspline2))>0) {
                    cout << "Error: k=" << k_bspline2 << " in "
                         << refinementbasis.mra.rangeI(j_bspline2) << " is missing. Supports: "
                         << refinementbasis.mra.phi.support(j_bspline2,k_bspline2) << " "
                         << refinementbasis.mra.phi.support(j_bspline1,k_bspline1) << endl;
                }
            }
            for (long k_bspline2=k_bspline2_last+1;
                      k_bspline2<=refinementbasis.mra.rangeI(j_bspline2).lastIndex(); ++k_bspline2) {
                if (overlap(refinementbasis.mra.phi.support(j_bspline1,k_bspline1),
                            refinementbasis.mra.phi.support(j_bspline2,k_bspline2))>0) {
                    cout << "Error: k=" << k_bspline2 << " in "
                         << refinementbasis.mra.rangeI(j_bspline2) << " is missing. Supports: "
                         << refinementbasis.mra.phi.support(j_bspline2,k_bspline2) << " "
                         << refinementbasis.mra.phi.support(j_bspline1,k_bspline1) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getWaveletNeighborsForBSpline(const PrimalBasis &basis, const RefinementBasis &refinementbasis)
{
    cout << " ******** Wavelet neighbors for BSpline **********" << endl;
    for (int refinement_j=refinementbasis.j0+4; refinement_j<=refinementbasis.j0+5; ++refinement_j) {
        for (int refinement_k =refinementbasis.mra.rangeI(refinement_j).firstIndex();
                 refinement_k<=refinementbasis.mra.rangeI(refinement_j).lastIndex(); ++refinement_k) {
            int j=0;
            long k_first=0L, k_last=0L;
            cout << "(" << refinement_j << "," << refinement_k << "): " << endl;
            refinementbasis.getWaveletNeighborsForBSpline(refinement_j, refinement_k, basis, j, k_first, k_last);
            cout << "(" << refinement_j << "," << refinement_k << "): "
                             << j << " , [" << k_first << "," << k_last << "], "
                             << basis.rangeJ(j) << endl;;

            for (long k=basis.rangeJ(j).firstIndex(); k<k_first; ++k) {
                if (overlap(refinementbasis.mra.phi.support(refinement_j,refinement_k),
                            basis.psi.support(j,k))>0) {
                    cout << "Error: k=" << k << " in " << basis.rangeJ(j) << " is missing."
                         << refinementbasis.mra.phi.support(refinement_j,refinement_k)
                         << " " << basis.psi.support(j,k) << endl;
                }
            }
            for (long k=k_last+1; k<basis.rangeJ(j).lastIndex(); ++k) {
                if (overlap(refinementbasis.mra.phi.support(refinement_j,refinement_k),
                            basis.psi.support(j,k))>0) {
                    cout << "Error: k=" << k << " in " << basis.rangeJ(j) << " is missing."
                         << refinementbasis.mra.phi.support(refinement_j,refinement_k)
                         << " " << basis.psi.support(j,k) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getScalingNeighborsForScaling(const PrimalBasis &basis)
{
    cout << " ******** Scaling neighbors for Scaling **********" << endl;
    for (int j_scaling1=basis.j0; j_scaling1<basis.j0+6; ++j_scaling1) {
        for (long k_scaling1= basis.mra.rangeI(j_scaling1).firstIndex();
                  k_scaling1<=basis.mra.rangeI(j_scaling1).lastIndex(); ++k_scaling1) {
            int j_scaling2=0;
            long k_scaling_first=0L, k_scaling_last=0L;
            basis.getScalingNeighborsForScaling(j_scaling1, k_scaling1, basis,
                                                j_scaling2, k_scaling_first, k_scaling_last);
            cout << "Scaling (" << j_scaling1 << "," << k_scaling1 << "): "
                                << j_scaling2 << " , [" << k_scaling_first << "," << k_scaling_last << "], "
                                << basis.mra.rangeI(j_scaling2) << endl;
            for (long k_scaling2=basis.mra.rangeI(j_scaling2).firstIndex();
                      k_scaling2<k_scaling_first; ++k_scaling2) {
                if (overlap(basis.mra.phi.support(j_scaling2,k_scaling2),
                            basis.mra.phi.support(j_scaling1,k_scaling1))>0) {
                    cout << "Error: k=" << k_scaling2 << " in " << basis.mra.rangeI(j_scaling2) << " is missing."
                         << basis.mra.phi.support(j_scaling2,k_scaling2)
                         << " " << basis.mra.phi.support(j_scaling1,k_scaling1) << endl;
                }
            }
            for (long k_scaling2=k_scaling_last+1;
                      k_scaling2<=basis.mra.rangeI(j_scaling2).lastIndex(); ++k_scaling2) {
                if (overlap(basis.mra.phi.support(j_scaling2,k_scaling2),
                            basis.mra.phi.support(j_scaling1,k_scaling1))>0) {
                    cout << "Error: k=" << k_scaling2 << " in " << basis.mra.rangeI(j_scaling2) << " is missing."
                         << basis.mra.phi.support(j_scaling2,k_scaling2)
                         << " " << basis.mra.phi.support(j_scaling1,k_scaling1) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getWaveletNeighborsForScaling(const PrimalBasis &basis)
{
    cout << " ******** Wavelet neighbors for Scaling **********" << endl;
    for (int j_scaling=basis.j0; j_scaling<basis.j0+6; ++j_scaling) {
        for (long k_scaling= basis.mra.rangeI(j_scaling).firstIndex();
                  k_scaling<=basis.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
            int j_wavelet=0;
            long k_wavelet_first=0L, k_wavelet_last=0L;
            basis.getWaveletNeighborsForScaling(j_scaling, k_scaling, basis,
                                                j_wavelet, k_wavelet_first, k_wavelet_last);
            cout << "Scaling (" << j_scaling << "," << k_scaling << "): "
                                << j_wavelet << " , [" << k_wavelet_first << "," << k_wavelet_last << "], "
                                << basis.rangeJ(j_wavelet) << endl;
            for (long k_wavelet=basis.rangeJ(j_wavelet).firstIndex();
                      k_wavelet<k_wavelet_first; ++k_wavelet) {
                if (overlap(basis.psi.support(j_wavelet,k_wavelet),
                            basis.mra.phi.support(j_scaling,k_scaling))>0) {
                    cout << "Error: k=" << k_wavelet << " in " << basis.rangeJ(j_wavelet) << " is missing."
                         << basis.psi.support(j_wavelet,k_wavelet)
                         << " " << basis.mra.phi.support(j_scaling,k_scaling) << endl;
                }
            }
            for (long k_wavelet=k_wavelet_last+1;
                      k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
                if (overlap(basis.psi.support(j_wavelet,k_wavelet),
                            basis.mra.phi.support(j_scaling,k_scaling))>0) {
                    cout << "Error: k=" << k_wavelet << " in " << basis.rangeJ(j_wavelet) << " is missing."
                         << basis.psi.support(j_wavelet,k_wavelet)
                         << " " << basis.mra.phi.support(j_scaling,k_scaling) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getBSplineNeighborsForWavelet(const PrimalBasis &basis, const RefinementBasis &refinementbasis)
{
    cout << " ******** BSpline neighbors for Wavelet **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<basis.j0+4; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
            int j_bspline=0;
            long k_bspline_first=0L, k_bspline_last=0L;
            basis.getBSplineNeighborsForWavelet(j_wavelet, k_wavelet, refinementbasis,
                                                j_bspline, k_bspline_first, k_bspline_last);
            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): "
                                << j_bspline << " , [" << k_bspline_first << "," << k_bspline_last << "], "
                                << refinementbasis.mra.rangeI(j_bspline) << endl;
            for (long k_bspline=refinementbasis.mra.rangeI(j_bspline).firstIndex();
                      k_bspline<k_bspline_first; ++k_bspline) {
                if (overlap(refinementbasis.mra.phi.support(j_bspline,k_bspline),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_bspline << " in " << refinementbasis.mra.rangeI(j_bspline) << " is missing."
                         << refinementbasis.mra.phi.support(j_bspline,k_bspline)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            for (long k_bspline=k_bspline_last+1;
                      k_bspline<=refinementbasis.mra.rangeI(j_bspline).lastIndex(); ++k_bspline) {
                if (overlap(refinementbasis.mra.phi.support(j_bspline,k_bspline),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_bspline << " in " << refinementbasis.mra.rangeI(j_bspline) << " is missing."
                         << refinementbasis.mra.phi.support(j_bspline,k_bspline)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getScalingNeighborsForWavelet(const PrimalBasis &basis)
{
    cout << " ******** Scaling neighbors for Wavelet **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<basis.j0+6; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
            int j_scaling=0;
            long k_scaling_first=0L, k_scaling_last=0L;
            basis.getScalingNeighborsForWavelet(j_wavelet, k_wavelet, basis,
                                                j_scaling, k_scaling_first, k_scaling_last);
            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): "
                                << j_scaling << " , [" << k_scaling_first << "," << k_scaling_last << "], "
                                << basis.mra.rangeI(j_scaling) << endl;
            for (long k_scaling=basis.mra.rangeI(j_scaling).firstIndex();
                      k_scaling<k_scaling_first; ++k_scaling) {
                if (overlap(basis.mra.phi.support(j_scaling,k_scaling),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_scaling << " in " << basis.mra.rangeI(j_scaling) << " is missing."
                         << basis.mra.phi.support(j_scaling,k_scaling)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            for (long k_scaling=k_scaling_last+1;
                      k_scaling<=basis.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
                if (overlap(basis.mra.phi.support(j_scaling,k_scaling),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_scaling << " in " << basis.mra.rangeI(j_scaling) << " is missing."
                         << basis.mra.phi.support(j_scaling,k_scaling)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getWaveletNeighborsForWavelet(const PrimalBasis &basis)
{
    cout << " ******** Wavelet neighbors for Wavelet **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<basis.j0+6; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
            int j_wavelet2=0;
            long k_wavelet_first=0L, k_wavelet_last=0L;
            basis.getWaveletNeighborsForWavelet(j_wavelet, k_wavelet, basis,
                                                j_wavelet2, k_wavelet_first, k_wavelet_last);
            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): "
                                << j_wavelet2 << " , [" << k_wavelet_first << "," << k_wavelet_last << "], "
                                << basis.rangeJ(j_wavelet2) << endl;
            for (long k_wavelet2=basis.rangeJ(j_wavelet2).firstIndex();
                      k_wavelet2<k_wavelet_first; ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet,k_wavelet),
                            basis.psi.support(j_wavelet2,k_wavelet2))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            for (long k_wavelet2=k_wavelet_last+1;
                      k_wavelet2<=basis.rangeJ(j_wavelet2).lastIndex(); ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet2,k_wavelet2),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getLowerWaveletNeighborsForWavelet(const PrimalBasis &basis)
{
    cout << " ******** Lower wavelet neighbors for Wavelet **********" << endl;
    for (int j_wavelet=basis.j0+1; j_wavelet<basis.j0+6; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
            int j_wavelet2=0;
            long k_wavelet_first=0L, k_wavelet_last=0L;
            basis.getLowerWaveletNeighborsForWavelet(j_wavelet, k_wavelet, basis,
                                                j_wavelet2, k_wavelet_first, k_wavelet_last);
            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): "
                                << j_wavelet2 << " , [" << k_wavelet_first << "," << k_wavelet_last << "], "
                                << basis.rangeJ(j_wavelet2) << endl;
            for (long k_wavelet2=basis.rangeJ(j_wavelet2).firstIndex();
                      k_wavelet2<k_wavelet_first; ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet,k_wavelet),
                            basis.psi.support(j_wavelet2,k_wavelet2))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            for (long k_wavelet2=k_wavelet_last+1;
                      k_wavelet2<=basis.rangeJ(j_wavelet2).lastIndex(); ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet2,k_wavelet2),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_getHigherWaveletNeighborsForWavelet(const PrimalBasis &basis)
{
    cout << " ******** Higher wavelet neighbors for Wavelet **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<basis.j0+6; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
            int j_wavelet2=0;
            long k_wavelet_first=0L, k_wavelet_last=0L;
            basis.getHigherWaveletNeighborsForWavelet(j_wavelet, k_wavelet, basis,
                                                j_wavelet2, k_wavelet_first, k_wavelet_last);
            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): "
                                << j_wavelet2 << " , [" << k_wavelet_first << "," << k_wavelet_last << "], "
                                << basis.rangeJ(j_wavelet2) << endl;
            for (long k_wavelet2=basis.rangeJ(j_wavelet2).firstIndex();
                      k_wavelet2<k_wavelet_first; ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet,k_wavelet),
                            basis.psi.support(j_wavelet2,k_wavelet2))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            for (long k_wavelet2=k_wavelet_last+1;
                      k_wavelet2<=basis.rangeJ(j_wavelet2).lastIndex(); ++k_wavelet2) {
                if (overlap(basis.psi.support(j_wavelet2,k_wavelet2),
                            basis.psi.support(j_wavelet,k_wavelet))>0) {
                    cout << "Error: k=" << k_wavelet2 << " in " << basis.rangeJ(j_wavelet2) << " is missing."
                         << basis.psi.support(j_wavelet2,k_wavelet2)
                         << " " << basis.psi.support(j_wavelet,k_wavelet) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void
test_precisionOfGenerators()
{
    long double OneDivSqrt2 = 0.707106781186547524401L;
    cout << 1./std::sqrt(2.L)-OneDivSqrt2 << " " << std::pow(2.L,-0.5L)-OneDivSqrt2
         << " " << pow2ih<long double>(-1)-OneDivSqrt2 << " " << pow2ih<double>(-1)-OneDivSqrt2 << endl;
    cout << 2.L/3.L << " " << _quadratic_refinement_left_evaluator0<long double>(2.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _quadratic_refinement_left_evaluator0<long double>(4.L/3.L,0) << endl << endl;
    cout << 4.L/9.L << " " << _cubic_refinement_left_evaluator0<long double>(1.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _cubic_refinement_right_evaluator0<long double>(4.L/3.L,0) << endl << endl;

}

/*
    MultiRefinementIntegral   refinement_integral(basis.refinementbasis,basis.refinementbasis);
    MultiWaveletIntegral multi_integral(basis,basis);
    DenseVectorLD *refCoeffs1, *refCoeffs2;
    T max_error = 0.L;
    for (int j1=0; j1<=J; ++j1) {
    for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
        for (int j2=0; j2<=J; ++j2) {
        for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
            int refinement_j1 = 0, refinement_j2 = 0;
            long refinement_k_first1 = 0L, refinement_k_first2 = 0L;
            refCoeffs1 = basis.psi.getRefinement(j1,k1,refinement_j1,refinement_k_first1);
            refCoeffs2 = basis.psi.getRefinement(j2,k2,refinement_j2,refinement_k_first2);
            T val = 0.L;
            for (int i1=(*refCoeffs1).firstIndex(); i1<=(*refCoeffs1).lastIndex(); ++i1) {
                for (int i2=(*refCoeffs2).firstIndex(); i2<=(*refCoeffs2).lastIndex(); ++i2) {
                    T integral_value = refinement_integral(refinement_j1,refinement_k_first1+i1,XBSpline,0,
                                                           refinement_j2,refinement_k_first2+i2,XBSpline,0);
                    val += (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value;
                }
            }
            if (j1==j2 && k1==k2) {
                cout << "(" << j1 << "," << k1 << "), (" << j2 << "," << k2 << "), Error: " << fabs(val-1.) << endl;
                max_error = std::max(max_error,fabs(val-1.));
            }
            else  {
                cout << "(" << j1 << "," << k1 << "), (" << j2 << "," << k2 << "), Error: " << fabs(val) << endl;
                max_error = std::max(max_error,fabs(val));
            }

        }
        }
    }
    }
    cout << "Max integration error: " << max_error << endl;
    */
