#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

/// FLENS typedefs: only required for the set up of right-hand side vectors
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Wavelet basis over an interval: Here we are exclusively (!!) considering $L_2$-orthonormal
///  multiwavelets.
typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

///  Definition of the (tensor product) wavelet preconditioner
typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying univariate bilinear form: As we are using $L_2$-orthonormal multiwavelets and
///  and are considering Poissons problem, we only need the univariate bilinear form
///  $a(v,w) = \int_0^1 v'(x) w'(x) dx$
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>      LaplaceOp1D;

///  Local operator in 1d for the above bilinear form
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D,LaplaceOp1D>          LocalOp1D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{A} \otimes \vec{\textrm{Id}}$
typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{\textrm{Id}} \otimes \vec{A}$
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

///  Aggregation of the two above two-dimensional operator in one class.
typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;

///  Righthandsides definitions: Here we are only considering right-hand side functions that can be
///  separated w.r.t.\ the coordinate direction, i.e., $f(x_1,x_2) = f_1(x_1) \otimes f_2(x)
///  or a sum of such functions.
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;

///  Multitree AWGM solver class.
typedef MultiTreeAWGM<Index2D,Basis2D,CompoundLocalOperator2D,
                      CompoundRhsIntegral2D,Preconditioner>         MultiTreeAWGM2D;

///  Some iterators we require for postprocessing.
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;


int example = 2;

T u1(T x)   {    return 1.; }

T u2(T y)   {    return 1.; }

T du1(T x)  {    return 0.; }

T du2(T y)  {    return 0.; }

T ddu1(T x) {    return -10.;   }

T ddu2(T y) {    return -10.;   }

long double EnergyErrorSquared = 0.L;

/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }

T u2(T y)   {    return y*y*(1-y)*(1-y); }

T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }

T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }

T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }

T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }

long double EnergyErrorSquared = 2.*(1.L/630.L * 2.L/105.L);
*/
T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index2D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }

    ///  Wavelet basis parameters
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    ///  Bulk chasing parameter
    T alpha = 0.7;

    ///  Relative tolerance for solving the finite-dimensional cg system in each iteration
    T gamma = 0.1;

    ///  Residual type: `Standard` refers to the construction proposed in Section 7.3
    const char* residualType = "standard";

    ///  Tree type we are using: `Sparse tree` refers to multitrees as introduced in Section 6.4.
    ///  `Graded tree` refers to the case where in one dimension, all overlappping wavelets need to
    ///  be included (see p.116)
    const char* treeType = "sparsetree"; //"gradedtree";

    ///  For this particular example, this is the only option. However, the multitree solver does
    ///  also work for other tensor product bases as it only requires an object that realizes
    ///  the matrix vector multiplication
    bool IsMW = true;

    ///  If required, ${\boldsymbol f}_\varepsilon$ can be computed and passed to the algorithms
    ///  This is only for testing the quality of the residual which can also be done with a forthcoming
    ///  program if we store (write) the coefficients to file
    bool compute_f_minus_Au_error = false;
    bool writeCoefficientsToFile = false;

    T eps   = 1e-5;
    Timer time;


    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Initialization of the operators
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,1.,1.,0.);

    /// Initialization of the right-hand side
    DenseVectorT sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    int order = 20;
    if (example==2) {  int order = 4+2*d; }
    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u2,sing_pts_y), fct_f2(f2,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral2D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data);
    AdaptiveSeparableRhsIntegral2D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data);
    CompoundRhsIntegral2D          F(rhs1,rhs2);

    Coefficients<Lexicographical,T,Index2D> f_eps(SIZEHASHINDEX2D);
    if (compute_f_minus_Au_error) {
        setUp_f_eps<T>(example, basis, Prec, f_eps,
                       rhs_u1_data, rhs_u1, rhs_u2_data, rhs_u2,
                       rhs_f1_data, rhs_f1, rhs_f2_data, rhs_f2);
    }

    ///  Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM2D multiTreeAWGM2D(basis2d, localOp2D, F, Prec/*, f_eps*/);
    multiTreeAWGM2D.setParameters(alpha, gamma, residualType, treeType, IsMW,
                                  /*compute_f_minus_Au_error,*/ writeCoefficientsToFile);

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    getSparseGridVector(basis2d,u,0,(T)0.2);

    ///  Set up strings for debugging respectively convergence measurement. In particular, storing
    ///  the coefficients in the indicated file allows testing the multitree based residual in
    ///  a forthcoming program.
    stringstream convfilename;
    convfilename << "conv_multitree_mw_awgm_poisson2d_" << example << "_" << argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_"
                 << treeType << ".dat";
    stringstream coefffilename;
    coefffilename << "coeff_multitree_mw_awgm_poisson2d_" << example << "_" << argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_" << treeType;

    ///  Calling the multitree AWGM solver
    multiTreeAWGM2D.cg_solve(u, eps, 100, 1e-2, EnergyErrorSquared, convfilename.str().c_str(),
                             coefffilename.str().c_str());

    ///  Plot the obtained solution, i.e., write point evaluations to a file for plotting with
    ///  gnuplot.
    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.1, "multiTreeAWGM_sol");

    return 0;
}

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index2D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2)
{
    int j0 = basis.j0;
    if (example==2) {
        for (int k=basis.mra.rangeIL(j0).firstIndex(); k<=basis.mra.rangeIL(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
        }
        for (int k=basis.mra.rangeIR(j0).firstIndex(); k<=basis.mra.rangeIR(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
        }

        for (int j=j0; j<=25; ++j) {
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            }
        }

        for (const_coeff1d_it it_x=rhs_u1_data.begin(); it_x!=rhs_u1_data.end(); ++it_x) {
            for (const_coeff1d_it it_y=rhs_u2_data.begin(); it_y!=rhs_u2_data.end(); ++it_y) {
                Index2D index((*it_x).first,(*it_y).first);
                f_eps[index] =  ( (*it_x).second * rhs_f2_data[(*it_y).first]
                               + rhs_f1_data[(*it_x).first] * (*it_y).second) * Prec(index);
            }
        }
    }
    std::cerr << "#Supp f_eps = " << f_eps.size() << std::endl;
}
