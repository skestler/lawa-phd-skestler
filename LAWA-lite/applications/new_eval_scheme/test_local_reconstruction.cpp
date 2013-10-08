/* TEST LOCAL RECONSTRUCTION
 *
 *  This examples calculates the local single reconstruction for wavelet interval bases. We consider
 *  Dijkema wavelet as well as L2-orthonormal multiwavelets.
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

/// We require special datastructures for univariate trees where coefficients are stored level-wise.
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

using namespace std;
using namespace lawa;

typedef double T;


///  Wavelet basis over an interval
//typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Typedef for iterators
typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;

///  Routine for constructing a random, nontrivial tree.
void
constructRandomTree(const PrimalBasis &basis, int J, TreeCoefficients1D<T> &LambdaTree);

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }

    srand (time(NULL));

    // wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);
    int j  = j0+5;
    int J  = atoi(argv[3]);

    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;

    /// The class LocalRefinement provides routines for the refinement of scaling functions
    /// $\phi_{j,k}$ and wavelets $\psi_{j,k}$ in terms of refinement B-Splines $\bar \phi_{j,k}$.
    LocalRefinement<PrimalBasis> LocalRefine(basis);

    /// Computing a vector of (multi-)scaling coefficients and a vector of (multi-)wavelet
    /// coefficients. Here, we consider a vector $v$ with
    /// $\mathrm{supp}\, v = \mathcal{I}_j \cup \mathcal{J}_j$, i.e., "full levels".
    CoefficientsByLevel<T> u_scaling1, u_wavelet1;
    for (int k=basis.mra.rangeI(j).firstIndex(); k<=basis.mra.rangeI(j).lastIndex(); ++k) {
        u_scaling1.map.operator[](k) = (T)rand()/RAND_MAX;
    }
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        u_wavelet1.map.operator[](k) = (T)rand()/RAND_MAX;
    }

    /// In case of multiscaling functions, we have to represent them first in a B-spline basis
    /// (refinement basis). Here, $\text{u_scaling1}$ contains the coefficients and $j$ indicates
    /// the level. The variables $\text{u_bspline1}$ and $\text{refinement_j_bspline}$ are called
    /// by reference and contain the corresponding coefficients in the refinement B-Spline basis on
    /// $\text{level refinement_j_bspline}$ so that
    /// $\sum_{v_k \in \text{u_scaling1}} v_k \phi_{j,k} =
    /// \sum_{\bar{d}_k \in \text{u_bspline1}} \bar{d}_k \bar{\phi}_{\text{refinement_j_bspline},k}$
    CoefficientsByLevel<T> u_bspline1;
    int refinement_j_bspline = 0;
    if (PrimalBasis::Cons==Multi && d>1) {
        LocalRefine.reconstructOnlyMultiScaling(u_scaling1, j, u_bspline1, refinement_j_bspline);
    }
    else {
        refinement_j_bspline = j;
        u_bspline1 = u_scaling1;
    }

    /// Now we compute the common representation of the refinement B-spline coefficient vector and
    /// the wavelet coefficient vector so that $\sum_{\bar{d}_k \in \text{u_bspline1}} \bar{d}_k \bar{\phi}_{\text{refinement_j_bspline},k}$
    /// $+ \sum_{c_k \in \text{u_wavelet1}} c_k \psi_{j,k} = \sum_{\bar{v}_k \in \text{u_loc_single_1}} \bar{v}_k \bar{\phi}_{\text{refinement_j},k}$
    CoefficientsByLevel<T> u_loc_single1;
    int refinement_j = 0;
    Timer time;
    time.start();
    LocalRefine.reconstruct(u_bspline1, refinement_j_bspline, u_wavelet1, j, u_loc_single1, refinement_j);
    time.stop();
    cout << "Local reconstruction took: " << time.elapsed() << endl;


    /// We validate our operations by comparing the error between our computed expansions. Here,
    /// $\text{val1}$ corresponds to the value of the multiscale representation, $\text{val2}$
    /// to the representation where the scaling functions is represented in terms of refinement
    /// B-Splines, and $\text{val3}$ corresponds to the representation by refinement B-Splines only.
    T max_error1 = 0., max_error2 = 0.;
    for (T x=0.; x<=1.; x+=pow2i<T>(-8-j)) {
        T val1 = 0.L, val2 = 0.L, val3 = 0.L;
        // no refinement
        for (const_coeffbylevel_it it=u_scaling1.map.begin(); it!=u_scaling1.map.end(); ++it) {
            val1 += (*it).second * basis.generator(XBSpline).operator()(x,j,(*it).first,0);
        }
        for (const_coeffbylevel_it it=u_wavelet1.map.begin(); it!=u_wavelet1.map.end(); ++it) {
            val1 += (*it).second * basis.generator(XWavelet).operator()(x,j,(*it).first,0);
        }
        // from multiscaling to bspline representation
        for (const_coeffbylevel_it it=u_bspline1.map.begin(); it!=u_bspline1.map.end(); ++it) {
            val2 += (*it).second * refinementbasis.generator(XBSpline).operator()(x,refinement_j_bspline,(*it).first,0);
        }
        for (const_coeffbylevel_it it=u_wavelet1.map.begin(); it!=u_wavelet1.map.end(); ++it) {
            val2 += (*it).second * basis.generator(XWavelet).operator()(x,j,(*it).first,0);
        }
        // single scale representation
        for (const_coeffbylevel_it it=u_loc_single1.map.begin(); it!=u_loc_single1.map.end(); ++it) {
            val3 += (*it).second * refinementbasis.generator(XBSpline).operator()(x,refinement_j,(*it).first,0);
        }
        max_error1 = std::max(max_error1,fabs(val1-val2));
        max_error2 = std::max(max_error2,fabs(val1-val3));
    }

    cout << "Reconstruction error for one step:" << max_error1 << " " << max_error2 << endl;


    /// Next, we intend to perform analogous operations as above on a tree. To this end, we
    /// construct a random tree on which we test the transformation to a local refinement B-Spline
    /// representation
    TreeCoefficients1D<T> u_tree(4096,basis.j0);
    Coefficients<Lexicographical,T,Index1D> u, u_loc_single;
    constructRandomTree(basis, J, u_tree);
    fromTreeCoefficientsToCoefficients(u_tree,u);
    IndexSet<Index1D> supp_u;
    supp_u = supp(u);
    for (const_set1d_it it=supp_u.begin(); it!=supp_u.end(); ++it) {
        u[*it] = (T)rand()/RAND_MAX;
    }

    /// The vector u contains the multilevel representation of a coefficient vector. We transform
    /// it to the local single scale representation.
    time.start();
    LocalRefine.reconstruct(u, j0, u_loc_single);
    time.stop();
    cout << "Local reconstruction took " << time.elapsed() << endl;

    /// We validate our operations by comparing the error between our computed expansions. Here,
    /// $\text{val1}$ corresponds to the value of the multiscale representation, $\text{val2}$
    /// to the representation by refinement B-Splines only.
    max_error1=0.;
    for (T x=0.; x<1.; x+=pow2i<T>(-6-J)) {
        T val1=0.L, val2=0.L;
        for (const_coeff1d_it it=u.begin(); it!=u.end(); ++it) {
            val1 += (*it).second *
                    basis.generator((*it).first.xtype).operator()(x,(*it).first.j,(*it).first.k,0);
        }
        for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
            val2 += (*it).second *
                    refinementbasis.mra.phi.operator()(x,(*it).first.j,(*it).first.k,0);
        }
        max_error1 = std::max(max_error1,fabs(val1-val2));
    }
    cout << "Reconstruction error for the whole tree: " << max_error1 << endl;

    /// Finally, we visualize the corresponding index sets. Visualizations of the corresponding
    /// index sets may then be created with Gnuplot.
    plotCoeff(u, basis, "coeff_multi_scale", false, true);
    plotCoeff(u_loc_single, basis, "coeff_local_single_scale", true, true);

    return 0;
}

void
constructRandomTree(const PrimalBasis &basis, int J, TreeCoefficients1D<T> &LambdaTree)
{
    /*
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaTree[0].map.operator[](k) = 0.;
    }
    for (int j=basis.j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            LambdaTree.bylevel[j-basis.j0+1].map.operator[](k) = 0.;
        }
    }
    */

    // First, we add some random values at random positions...
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaTree[0].map.operator[](k) = 0.;
    }
    for (int j=basis.j0; j<=J; ++j) {
        int random_k1 = rand() % basis.cardJ(j) + 1;
        LambdaTree.bylevel[j-basis.j0+1].map.operator[](random_k1) = 0.;
        int random_k2 = rand() % basis.cardJ(j) + 1;
        LambdaTree.bylevel[j-basis.j0+1].map.operator[](random_k2) = 0.;
    }

    // ... and then, going from top to bottom, we add indices to make the index set a tree.
    for (int i=J-basis.j0+1; i>=2; --i) {
        CoefficientsByLevel<T> *currentlevel;
        currentlevel = &(LambdaTree.bylevel[i]);
        int j=basis.j0+i-1;
        for (const_coeffbylevel_it it=(*currentlevel).map.begin(); it!=(*currentlevel).map.end(); ++it) {
            long k = (*it).first;
            long k_first = basis.rangeJ(j-1).firstIndex();
            long k_last  = basis.rangeJ(j-1).lastIndex();
            for (int k1=k_first; k1<=k_last; ++k1) {
                if (overlap(basis.psi.support(j,k),basis.psi.support(j-1,k1))>0) {
                    LambdaTree.bylevel[i-1].map.operator[](k1) = 0.;
                }
            }
        }
    }
}

/*
 * IntegralRefinentBasis integral_refinement(refinementbasis,refinementbasis);
    int N = refinementbasis.mra.cardI(j);
    DenseMatrixT A(N,N);
    int offset = refinementbasis.mra.rangeI(j).firstIndex()-1;
    for (int k_row=refinementbasis.mra.rangeI(j).firstIndex(); k_row<=refinementbasis.mra.rangeI(j).lastIndex(); ++k_row) {
        for (int k_col=refinementbasis.mra.rangeI(j).firstIndex(); k_col<=refinementbasis.mra.rangeI(j).lastIndex(); ++k_col) {
            A(k_row-offset,k_col-offset) = integral_refinement(j,k_row,XBSpline,1, j,k_col,XBSpline,1);
        }
    }
    DenseVector<Array<T> > wr(N), wi(N);
    DenseMatrixT vl,vr;
    ev(false, false, A, wr, wi, vl, vr);
    T cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
    cout << "Eigenvalues for A_" << j <<", kappa = " << CB/cB << ", cA = " << cB << ", CA = " << CB << endl;
 */
