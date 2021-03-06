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
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis> Basis3D;

///  Definition of the (tensor product) wavelet preconditioner
typedef OptimizedH1Preconditioner3D<T,Basis3D>                      Preconditioner3D;

///  Underlying univariate bilinear form: As we are using $L_2$-orthonormal multiwavelets and
///  and are considering Poissons problem, we only need the univariate bilinear form
///  $a(v,w) = \int_0^1 v'(x) w'(x) dx$
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>      LaplaceOp1D;

///  Local operator in 1d for the above bilinear form
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D, LaplaceOp1D>         LocalOp1D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{A} \otimes \vec{\textrm{Id}} \otimes \vec{\textrm{Id}}$
typedef UniDirectionalLocalOperator<Index3D,XOne,LocalOp1D,
                                            NotXOne,Index2D>        UniDirectionalLocalOpXOne3D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{\textrm{Id}} \otimes \vec{A} \otimes \vec{\textrm{Id}}$
typedef UniDirectionalLocalOperator<Index3D,XTwo,LocalOp1D,
                                            NotXTwo,Index2D>        UniDirectionalLocalOpXTwo3D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{\textrm{Id}} \otimes \vec{\textrm{Id}} \otimes \vec{A}$
typedef UniDirectionalLocalOperator<Index3D,XThree,LocalOp1D,
                                            NotXThree,Index2D>      UniDirectionalLocalOpXThree3D;

///  Aggregation of the two above two-dimensional operator in one class.
typedef CompoundLocalOperator<Index3D,
                              UniDirectionalLocalOpXOne3D,
                              UniDirectionalLocalOpXTwo3D,
                              UniDirectionalLocalOpXThree3D>        CompoundLocalOperator3D;

///  Righthandsides definitions: Here we are only considering right-hand side functions that can be
///  separated w.r.t.\ the coordinate direction, i.e., $f(x_1,x_2) = f_1(x_1) \otimes f_2(x)
///  or a sum of such functions.
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index3D,Rhs1D,Rhs1D,Rhs1D>           AdaptiveSeparableRhsIntegral3D;
typedef CompoundRhs<T,Index3D,
                    AdaptiveSeparableRhsIntegral3D,
                    AdaptiveSeparableRhsIntegral3D,
                    AdaptiveSeparableRhsIntegral3D>                 CompoundRhsIntegral3D;

///  Multitree AWGM solver class.
typedef MultiTreeAWGM<Index3D,Basis3D,CompoundLocalOperator3D,
                      CompoundRhsIntegral3D,Preconditioner3D>       MultiTreeAWGM3D;

///  Some iterators we require for postprocessing.
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef IndexSet<Index3D>::const_iterator                           const_set3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator           coeff3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator     const_coeff3d_it;

int example = 2;
T u1(T x)   {    return 1.; }
T u2(T y)   {    return 1.; }
T u3(T z)   {    return 1.; }
T du1(T x)  {    return 0.; }
T du2(T y)  {    return 0.; }
T du3(T z)  {    return 0.; }
T ddu1(T x) {    return -33-(T)(1.L/3.L); }
T ddu2(T y) {    return -33-(T)(1.L/3.L); }
T ddu3(T z) {    return -33-(T)(1.L/3.L); }

long double EnergyErrorSquared = 14.20158453089639L*14.20158453089639L;

/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }
T u2(T y)   {    return y*y*(1-y)*(1-y); }
T u3(T z)   {    return z*z*(1-z)*(1-z); }
T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }
T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }
T du3(T z)  {    return 2*z*(1-z)*(1-z)-2*z*z*(1-z); }
T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }
T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }
T ddu3(T z) {    return 2*(1-z)*(1-z) - 8*z*(1-z) + 2*z*z; }

long double EnergyErrorSquared = 3.*(1.L/630.L * 1.L/630.L * 2.L/105.L);
*/
T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T f3(T z)   {   return -ddu3(z); }

T sol(T x, T y, T z) {   return u1(x) * u2(y) * u3(z); }

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
    const char* treeType = "sparsetree";//"gradedtree";

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
    Basis3D basis3d(basis,basis,basis);

    /// Initialization of the operators
    LaplaceOp1D                     laplaceOp1D(basis);
    LocalOp1D                       localOp1D(basis,basis,refinementbasis.LaplaceOp1D, laplaceOp1D);
    UniDirectionalLocalOpXOne3D     uniDirectionalOpXOne3D(localOp1D);
    //uniDirectionalOpXOne3D.setParameters(J, 49157, 6151);
    UniDirectionalLocalOpXTwo3D     uniDirectionalOpXTwo3D(localOp1D);
    //uniDirectionalOpXTwo3D.setParameters(J, 49157, 6151);
    UniDirectionalLocalOpXThree3D   uniDirectionalOpXThree3D(localOp1D);
    //uniDirectionalOpXThree3D.setParameters(J, 49157, 6151);
    CompoundLocalOperator3D         localOp3D(uniDirectionalOpXOne3D,uniDirectionalOpXTwo3D,
                                              uniDirectionalOpXThree3D);

    /// Initialization of preconditioner
    Preconditioner3D  Prec(basis3d,1.,1.,1.,0.);

    /// Initialization of the right-hand side
    DenseVectorT sing_pts_x, sing_pts_y, sing_pts_z;
    DenseMatrixT no_deltas, deltas_x, deltas_y, deltas_z;
    int order = 20;
    if (example == 2) order = 2*d;
    if (example == 3) order = 4+2*d;
    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u1,sing_pts_y), fct_f2(f2,sing_pts_y);
    Function<T>                    fct_u3(u3,sing_pts_z), fct_f3(f3,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u3(basis, fct_u3, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f3(basis, fct_f3, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D),
                                            rhs_u3_data(SIZEHASHINDEX1D),
                                            rhs_f3_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral3D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs3(rhs_u1, rhs_u1_data, rhs_u2, rhs_u2_data,
                                        rhs_f3, rhs_f3_data);
    CompoundRhsIntegral3D          F(rhs1,rhs2,rhs3);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM3D multiTreeAWGM3D(basis3d, localOp3D, F, Prec);
    multiTreeAWGM3D.setParameters(alpha, gamma, residualType, treeType, IsMW, writeCoefficientsToFile);

    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    getSparseGridVector(basis3d,u,0,(T)0.2);

    ///  Set up strings for debugging respectively convergence measurement. In particular, storing
    ///  the coefficients in the indicated file allows testing the multitree based residual in
    ///  a forthcoming program.
    stringstream convfilename;
    convfilename << "conv_multitree_mw_awgm_poisson3d_" << example << "_" <<argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_"
                 << treeType << ".dat";
    stringstream coefffilename;
    coefffilename << "coeff_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_" << treeType;

    ///  Calling the multitree AWGM solver
    multiTreeAWGM3D.cg_solve(u, eps, 100, 1e-2, EnergyErrorSquared,
                             convfilename.str().c_str(), coefffilename.str().c_str());

    return 0;
}


