/*
 * POISSON PROBLEM 1D
 *
 *  This example calculates a poisson problem with constant forcing f on the 
 *  one-dimensional domain [0,1], i.e.
 *          - u'' = f on (0,1) , u(0) = u(1) = 0.
 *  The solution is obtained using a uniform Wavelet-Galerkin method with a
 *  diagonal scaling preconditioner.
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;


typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef HelmholtzOperator1D<T, PrimalBasis>                         HelmholtzOp;
typedef DiagonalMatrixPreconditioner1D<T, PrimalBasis, HelmholtzOp> DiagonalPrec;
typedef RHSWithPeaks1D<T, PrimalBasis>                              Rhs;

T
rhs_f(T /*x*/)
{
    return 1.;
}

void
printU(const DenseVectorT u, const PrimalBasis& basis, const int J, 
       const char* filename, const double deltaX=1./128.)
{
    ofstream file(filename);
    for(double x = 0; x <= 1.; x += deltaX){
        file << x << " " << evaluate(basis,J, u, x, 0) << endl; 
    }
    file.close();
}

int main()
{
    int d = 2;      // (d,d_)-wavelets
    int d_ = 2;
    int j0 = 2;     // minimal level
    int J = 5;      // maximal level
    
    PrimalBasis basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    
    HelmholtzOp  a(basis, 0);
    DiagonalPrec p(a);
    
    DenseVectorT singPts;                      // singular points of the rhs forcing function: here none
    FullColMatrixT deltas;                     // peaks (and corresponding scaling coefficients): here none
    Function<T> F(rhs_f, singPts);             // Function object (wraps a function and its singular points)
    Rhs rhs(basis, F, deltas, 4, false, true); // RHS: specify integration order for Gauss quadrature (here: 4) and
                                               //      if there are singular parts (false) and/or smooth parts (true) 
                                               //      in the integral
    
    Assembler1D<T, PrimalBasis> assembler(basis);
    SparseMatrixT   A = assembler.assembleStiffnessMatrix(a, J);
    DiagonalMatrixT P = assembler.assemblePreconditioner(p, J);
    DenseVectorT    f = assembler.assembleRHS(rhs, J);
    
    DenseVectorT u(basis.mra.rangeI(J));
    
    cout << pcg(P, A, u, f) << " pcg iterations" << endl;
    
    printU(u, basis, J, "u.txt");
    
    return 0;
}
