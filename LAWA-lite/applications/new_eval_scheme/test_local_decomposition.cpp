/* TEST LOCAL DECOMPOSITION
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

///  Typedefs for FLENS and precision. Important: `long double` is only available for $L_2$-
///  orthonormal constructions!
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Wavelet basis over an interval and corresponding refinement B-Spline basis
//typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Wavelet integrals w.r.t. wavelet basis and refinement B-Splines
typedef IntegralF<Gauss,PrimalBasis>                                IntegralFBasis;
typedef IntegralF<Gauss,RefinementBasis>                            IntegralFRefinementBasis;

typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;

/// The function $f$ we use for this test.
T
f(T x) { return exp(x); }

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=3) {
        cout << "Usage: " << argv[0] << " d j0" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);

    /// Initialization of the level $j$ in our documentation.
    int j_wavelet  = j0+5;
    int j_scaling  = j_wavelet;

    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, j0);           // For L2_orthonormal MW bases
    PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;

    /// Integral initialization
    DenseVectorT                singPts;
    Function<T>                 f_fct(f,singPts);
    IntegralFBasis              integralf_basis(f_fct,basis);
    IntegralFRefinementBasis    integralf_refinementbasis(f_fct,refinementbasis);
    integralf_basis.quadrature.setOrder(30);
    integralf_refinementbasis.quadrature.setOrder(30);

    /// Local refinement initialization
    LocalRefinement<PrimalBasis> LocalRefine(basis);


    /// Get $j + \bar{j}+1$.
    int j_refinement = basis.psi.getRefinementLevel(j_wavelet);

    /// Computing a vector of refinement B-spline coefficients, i.e., $\mathbf{f}$.
    CoefficientsByLevel<T> f_loc_single;
    for (int k= refinementbasis.mra.rangeI(j_refinement).firstIndex();
             k<=refinementbasis.mra.rangeI(j_refinement).lastIndex(); ++k) {
        f_loc_single.map.operator[](k) = integralf_refinementbasis(j_refinement,k,XBSpline,0);
    }

    /// Initializing output vectors for (multi-)wavelets and (multi-)scaling functions. This is
    /// required as only existing entries in the corresponding hash maps are processed!
    CoefficientsByLevel<T> f_bspline, f_wavelet, f_scaling;
    for (int k =basis.rangeJ(j_wavelet).firstIndex(); k<=basis.rangeJ(j_wavelet).lastIndex(); ++k) {
        f_wavelet.map.operator[](k) = 0.;
    }
    for (int k =refinementbasis.mra.rangeI(j_refinement-1).firstIndex();
             k<=refinementbasis.mra.rangeI(j_refinement-1).lastIndex(); ++k) {
        f_bspline.map.operator[](k) = 0.;
    }
    if (PrimalBasis::Cons==Multi) {
        for (int k =basis.mra.rangeI(j_scaling).firstIndex();
                 k<=basis.mra.rangeI(j_scaling).lastIndex(); ++k) {
            f_scaling.map.operator[](k) = 0.;
        }
    }

    /// Decompose the refinement B-spline coefficient vector, i.e., compute $\bar{\mathbf{g}}$ and
    /// $\mathbf{h}$.
    LocalRefine.decompose_(f_loc_single, f_bspline, j_refinement-1, f_wavelet, j_wavelet);

    cout << "************** Wavelet Decompositions **************" << endl;
    /// Check if wavelet coefficients in $\mathbf{h}$ are correct.
    for (const_coeffbylevel_it it=f_wavelet.map.begin(); it!=f_wavelet.map.end(); ++it) {
        T val1 = integralf_basis(j_wavelet,(*it).first,XWavelet,0);
        T val2 = (*it).second;
        cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
             << fabs(val1-val2)/fabs(val1) << endl;
    }
    cout << endl;

    cout << "************** B-spline Decompositions *************" << endl;
    /// Check if refinement B-Spline coefficients in $\bar{\mathbf{g}}$ are correct.
    for (const_coeffbylevel_it it=f_bspline.map.begin(); it!=f_bspline.map.end(); ++it) {
        T val1 = integralf_refinementbasis(j_refinement-1,(*it).first,XBSpline,0);
        T val2 = (*it).second;
        cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
             << fabs(val1-val2)/fabs(val1) << endl;
    }
    cout << endl;

    /// In case of multiscaling functions, we have to further decompose the refinement B-splines,
    /// i.e., compute $\mathbf{g}$.
    if (PrimalBasis::Cons==Multi) {
        LocalRefine.decompose_OnlyMultiScaling(f_bspline, f_scaling, j_scaling);

        /// Check if refinement scaling function coefficients in $\mathbf{g}$ are correct.
        cout << "************** Scaling Decompositions *************" << endl;
        for (const_coeffbylevel_it it=f_scaling.map.begin(); it!=f_scaling.map.end(); ++it) {
            T val1 = integralf_basis(j_scaling,(*it).first,XBSpline,0);
            T val2 = (*it).second;
            cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
                 << fabs(val1-val2)/fabs(val1) << endl;
        }
    }


    return 0;
}

