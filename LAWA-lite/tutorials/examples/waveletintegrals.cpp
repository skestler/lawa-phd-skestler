/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Typedef for double precision
typedef double T;

///  Typedefs for Flens data types:
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

///  Typedefs for basis constructions (here for `Dijkema` and `Orthogonal` basis constructions):
///     Dijkema Basis over an interval
typedef Basis<T,Primal,Interval,Dijkema> Basis1D;
///     L2 orthonormal Basis over an interval
//typedef Basis<T,Orthogonal,Interval,Multi> Basis1D;

///  Typedefs for integral:
///     Integral $\int \psi_{j_1,k_1}(x) \psi_{j_2,k_2}(x) dx$ to be computed by `Gauss-Legendre`
///     quadrature
typedef Integral<Gauss,Basis1D,Basis1D> Integral1D;
///     Integral $\int f(x) \psi_{j_2,k_2}(x) dx$ to be computed by `Gauss-Legendre`
///     quadrature
typedef IntegralF<Gauss,Basis1D> IntegralF1D;
///     Integral $\int a(x) \psi_{j_1,k_1}(x) \psi_{j_2,k_2}(x) dx$ to be computed by
///     `Gauss-Legendre` quadrature
typedef IntegralF<Gauss,Basis1D,Basis1D> IntegralFF1D;

T
f(T x) {
    return exp(x);
}

T
a(T x) {
    if (x<1./3.)    return  1.;
    else            return -1.;
}

int main (int argc, char *argv[]) {

    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(1);
    }

    /// Wavelet basis parameters:
    int d  = atoi(argv[1]);     // polynomial order
    int d_ = atoi(argv[2]);     // vanishing moments
    int j0 = atoi(argv[3]);     // minimal level
    int J  = atoi(argv[4]);     // maximum level

    Basis1D basis(d,d_,j0);
    //Basis1D basis(d,j0);
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();

    /// Creating an instance of a function object
    ///     Store the singular values of function of interest in FLENS vector
    DenseVectorT singularPoints_f;
    DenseVectorT singularPoints_a(1);
    singularPoints_a = 1./3.;
    ///     Corresponding function objects
    Function<T> function_f(f, singularPoints_f);
    Function<T> function_a(a, singularPoints_a);

    /// Wavelet integrals
    Integral1D   integral(basis,basis);
    IntegralF1D  integralF(function_f,basis);
    IntegralFF1D integralFF(function_a,basis,basis);

    /// Order of integration
    ///     We set the order of the Gauss-Legendre quadrature to 10
    integralF.quadrature.setOrder(10);
    integralFF.quadrature.setOrder(10);

    /// Calling the integral routines
    int   j1 = j0 + 3, j2 = j0 + 5;
    long  k1 = 4, k2 = 6;
    XType xtype1 = XWavelet, xtype2 = XWavelet;
    int   deriv1 = 0, deriv2 = 0;

    cout << "wavelet vs. wavelet: " << integral(j1,k1,xtype1,deriv1, j2,k2,xtype2,deriv2)   << endl;
    cout << "function vs. wavelet: " << integralF(j1,k1,xtype1,deriv1, j2,k2,xtype2,deriv2) << endl;
    cout << "function vs. wavelet vs. wavelet: " << integralFF(j1,k1,xtype1,deriv1,
                                                               j2,k2,xtype2,deriv2)  << endl;

    return 0;
}

