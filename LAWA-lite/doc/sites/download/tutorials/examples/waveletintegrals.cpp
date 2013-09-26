#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

typedef Basis<T,Primal,Interval,Dijkema> Basis1D;
//typedef Basis<T,Orthogonal,Interval,Multi> Basis1D;

typedef Integral<Gauss,Basis1D,Basis1D> Integral1D;
typedef IntegralF<Gauss,Basis1D> IntegralF1D;
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

    int d  = atoi(argv[1]);     // polynomial order
    int d_ = atoi(argv[2]);     // vanishing moments
    int j0 = atoi(argv[3]);     // minimal level
    int J  = atoi(argv[4]);     // maximum level

    Basis1D basis(d,d_,j0);
    //Basis1D basis(d,j0);
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();

    DenseVectorT singularPoints_f;
    DenseVectorT singularPoints_a(1);
    singularPoints_a = 1./3.;
    Function<T> function_f(f, singularPoints_f);
    Function<T> function_a(a, singularPoints_a);

    Integral1D   integral(basis,basis);
    IntegralF1D  integralF(function_f,basis);
    IntegralFF1D integralFF(function_a,basis,basis);

    integralF.quadrature.setOrder(10);
    integralFF.quadrature.setOrder(10);

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
