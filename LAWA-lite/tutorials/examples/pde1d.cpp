/* POISSON PROBLEM 1D
 *
 *  This example calculates a poisson problem with constant forcing f on the
 *  one-dimensional domain [0,1], i.e.
 *          - u'' = f on (0,1) , u(0) = u(1) = 0.
 *  The solution is obtained using a uniform Wavelet-Galerkin method with a
 *  diagonal scaling preconditioner.
 */

/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:
///     Primal Basis over an interval, using Dijkema construction
//typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
//typedef Basis<T, Primal, Interval, SparseMulti>                         PrimalBasis;
typedef Basis<T, Orthogonal, Interval, Multi>                         PrimalBasis;

///     PDE-Operator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x) \int(b(x) v \cdot u') + \int(a(x) v \cdot u)$
typedef WeightedPDEOperator1D<T, PrimalBasis>                       PDEOp;

///     Preconditioner: diagonal scaling with norm of operator
typedef DiagonalMatrixPreconditioner1D<T, PrimalBasis, PDEOp>       Prec;

///     Right Hand Side (RHS): basic 1D class for rhs integrals of the form $\int(f \cdot v)$,
///     possibly with additional peak contributions (not needed here)
typedef RHSWithPeaks1D<T, PrimalBasis>                              Rhs;

/// Reference solution
T
u_f(T x)
{
    //return x*(x-1);
    return x*(x-1)*std::exp(x);
}

T
u_x_f(T x)
{
    //return 2*x-1;
    return (x*x+x-1.)*std::exp(x);
}

T
u_xx_f(T x)
{
    //return 2;
    return (x*x+3*x)*std::exp(x);
}

/// Non-constant reaction term
T
a_f(T x)
{
    return 1.;
    //return 1+exp(x);
}

/// Non-constant convection term
T
b_f(T x)
{
    //return 1+exp(x);
    return 1.;
}

/// constant diffusion term
T
c_f(T x)
{
    return 1.;  // For this example, this term needs to be constant in order to get a simple
                // representation for the rhs (see below).
}

/// Forcing function of the form `T f(T x)` - here a constant function
T
rhs_f(T x)
{
    return -c_f(x)*u_xx_f(x) + b_f(x)*u_x_f(x) + a_f(x)*u_f(x);
}

/// Auxiliary function to print solution values, generates `.txt`-file with
/// columns: `x u(x)`
void
printU(const DenseVectorT u, const PrimalBasis& basis, const int J,
       const char* filename, const double deltaX=1./128.)
{
    ofstream file(filename);
    for(double x = 0; x <= 1.; x += deltaX){
        file << x << " " << u_f(x) << " " << u_x_f(x)
                   << " " << evaluate(basis,J, u, x, 0) << " " << evaluate(basis,J, u, x, 1) << endl;
    }
    file.close();
}

/// Auxiliary function to print solution values, generates `.txt`-file with
/// columns: `x u(x)`
void
H1errorU(const DenseVectorT u, const PrimalBasis& basis, const int j, long double &L2error, long double &H1error,
         const double deltaX=1./128.)
{
    L2error = 0.;
    H1error = 0.;
    long double H1seminormerror = 0.;
    long double diff_u     = u_f(0.)     - evaluate(basis,j, u, 0., 0);
    long double diff_u_x   = u_x_f(0.)   - evaluate(basis,j, u, 0., 1);
    L2error         += 0.5*diff_u*diff_u;
    H1seminormerror += 0.5*diff_u_x*diff_u_x;
    for(double x = deltaX; x <= 1.-deltaX; x += deltaX) {
        diff_u   = u_f(x)   - evaluate(basis,j, u, x, 0);
        diff_u_x = u_x_f(x) - evaluate(basis,j, u, x, 1);
        L2error += diff_u*diff_u;
        H1seminormerror += diff_u_x*diff_u_x;
    }
    diff_u     = u_f(1.)     - evaluate(basis,j, u, 1., 0);
    diff_u_x   = u_x_f(1.)   - evaluate(basis,j, u, 1., 1);
    L2error         += 0.5*diff_u*diff_u;
    H1seminormerror += 0.5*diff_u_x*diff_u_x;
    H1error = L2error + H1seminormerror;
    L2error*= deltaX;
    H1error*= deltaX;
    L2error = std::sqrt(L2error);
    H1error = std::sqrt(H1error);
}

int main(int argc, char*argv[])
{
    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    /// wavelet basis parameters:
    int d = atoi(argv[1]);
    int d_ =atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J = atoi(argv[4]);

    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, d_, j0); // For biorthogonal wavelet bases
    PrimalBasis basis(d, j0);       // For L2_orthonormal and special MW bases
    basis.enforceBoundaryCondition<DirichletBC>();

    /// Operator initialization
    int order = 20; // quadrature order
    DenseVectorT a_singPts, b_singPts, c_singPts;
    Function<T> a(a_f, a_singPts);
    Function<T> b(b_f, b_singPts);
    Function<T> c(c_f, c_singPts);
    PDEOp       op(basis, a, b, c, order);
    Prec        p(op);

    /// Righthandside initialization
    DenseVectorT singPts;                      // singular points of the rhs forcing function: here none
    DenseMatrixT deltas;                     // peaks (and corresponding scaling coefficients): here none
    Function<T> F(rhs_f, singPts);             // Function object (wraps a function and its singular points)
    Rhs rhs(basis, F, deltas, 40, false, true); // RHS: specify integration order for Gauss quadrature (here: 4) and
                                               //      if there are singular parts (false) and/or smooth parts (true)
                                               //      in the integral

    ofstream file("pde1d_convergence.txt");
    for (int j=j0; j<=J; ++j) {
        /// Assembler: assemble the problem components
        Assembler1D<T, PrimalBasis> assembler(basis);
        SparseMatrixT   A = assembler.assembleStiffnessMatrix(op, j);
        DiagonalMatrixT P = assembler.assemblePreconditioner(p, j);
        DenseVectorT    f = assembler.assembleRHS(rhs, j);

        DenseMatrixT A_dense;
        densify(cxxblas::NoTrans,A,A_dense);
        int N = A_dense.numRows();

        for (int i=1; i<=N; ++i) {
            for (int j=1; j<=N; ++j) {
                A_dense(i,j) *= P._diag(i)*P._diag(j);
            }
        }
        DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());

        DenseVector<Array<T> > wr(N), wi(N);
        DenseMatrixT vl,vr;
        ev(false, false, A_dense, wr, wi, vl, vr);
        T cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
        for (int i=1; i<=wr.lastIndex(); ++i) {
            cB = std::min(cB,wr(i));
            CB = std::max(CB,wr(i));
        }
        cout << "Eigenvalues for A_" << j <<", kappa = " << CB/cB
             << ", cA = " << cB << ", CA = " << CB << endl;

        spy(A,"A");

        //getchar();

        /// Initialize empty solution vector
        DenseVectorT u(basis.mra.rangeI(j));

        /// Solve problem using pgmres
        cout << pgmres(P, A, u, f, 1e-12) << " pgmres iterations" << endl;

        /// Compute errors
        long double L2error=0.L, H1error=0.L;
        H1errorU(u, basis, j, L2error, H1error, 0.0000001/*pow2i<T>(-j-15)*/);
        file << basis.mra.cardI(j) << " " << L2error << " " << H1error << endl;
        cout << basis.mra.cardI(j) << " " << L2error << " " << H1error << endl;

        /// Print solution to file "u.txt"
        printU(u, basis, j, "u.txt", pow2i<T>(-j-10));
        std::cerr << "Loop finished." << std::endl;
    }
    file.close();
    return 0;
}
