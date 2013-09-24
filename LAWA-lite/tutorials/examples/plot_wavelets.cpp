/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Typedef for double precision
typedef double T;

///  Typedefs for basis construction:
///     Dijkema Basis over an interval
typedef Basis<T,Primal,Interval,Dijkema> Basis1D;
///     L2 orthonormal Basis over an interval
//typedef Basis<T,Orthogonal,Interval,Multi> Basis1D;

//typedef Basis<T,Orthogonal,R,Multi> Basis1D;
//typedef Basis<T,Orthogonal,R,MultiRefinement> Basis1D;

///  Typedefs for integral:
typedef Integral<Gauss,Basis1D,Basis1D> IntegralBasis1DvsBasis1D;

template<Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Interval,Cons>& basis, int j, XType type);

template<Construction Cons>
Range<int>
getRange(const Basis<T,Orthogonal,Interval,Cons>& basis, int j, XType type);

template<Construction Cons>
Range<int>
getRange(const Basis<T,Orthogonal,R,Cons>& basis, int j, XType type);

template<Construction Cons>
Range<int>
getRange(const Basis<T,Dual,Interval,Cons>& basis, int j, XType type);


int main (int argc, char *argv[]) {

    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    cout.precision(20);
    /// wavelet basis parameters:
    int d = atoi(argv[1]);
    int d_ =atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J = atoi(argv[4]);

    Basis1D basis(d,d_,j0);
    //Basis1D basis(d,j0);
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();
    //T a = -5., b = 5.;
    T a = 0., b = 1.;


    IntegralBasis1DvsBasis1D integral(basis,basis);

    Range<int> scalingrange = getRange(basis, j0, XBSpline);
    /// Plot scaling function and wavelets

    for (int k=scalingrange.firstIndex(); k<=scalingrange.lastIndex(); ++k) {
        stringstream filename;
        filename << "scaling_" << k << ".txt";
        ofstream plotfile_scaling(filename.str().c_str());
        for (T x=basis.mra.phi.support(j0,k).l1; x<=basis.mra.phi.support(j0,k).l2; x+=pow2i<T>(-6)) {
            plotfile_scaling << x << " " << basis.generator(XBSpline)(x,j0,k,0) << endl;
            if (x==0.) {
                cout << "supp phi_{"<< j0 << "," << k << "} = " << basis.generator(XBSpline).support(j0,k) << endl;
                cout << "||phi_{" << j0 << "," << k << "}||^2_L2 = "
                     << basis.generator(XBSpline).getL2Norm(j0,k) << " "
                     << std::sqrt(integral(j0,k,XBSpline,0, j0,k,XBSpline,0)) << endl;
                cout << "||d/dx  phi_{" << j0 << "," << k << "}||^2_L2 = "
                     << basis.generator(XBSpline).getH1SemiNorm(j0,k) << " "
                     << std::sqrt(integral(j0,k,XBSpline,1, j0,k,XBSpline,1)) << endl << endl;
            }
        }
    }

    for (int j=j0; j<=J-1; ++j) {
        Range<int> waveletrange = getRange(basis, j, XWavelet);
        for (int k=waveletrange.firstIndex(); k<=waveletrange.lastIndex(); ++k) {
            stringstream filename;
            filename << "wavelet_" << j << "_" << k << ".txt";
            ofstream plotfile_wavelet(filename.str().c_str());
            for (T x=basis.psi.support(j,k).l1; x<=basis.psi.support(j,k).l2; x+=pow2i<T>(-6)) {
                plotfile_wavelet << x << " " << basis.generator(XWavelet)(x,j,k,0) << endl;
                if (x==0.) {
                    cout << "supp psi_{"<< j << "," << k << "} = " << basis.generator(XWavelet).support(j,k) << endl;
                    cout << "||psi_{" << j << "," << k << "}||^2_L2 = "
                         << basis.generator(XWavelet).getL2Norm(j,k) << " "
                         << std::sqrt(integral(j,k,XWavelet,0, j,k,XWavelet,0)) << endl;
                    cout << "||d/dx  psi_{" << j << "," << k << "}||^2_L2 = "
                         << basis.generator(XWavelet).getH1SemiNorm(j,k) << " "
                         << std::sqrt(integral(j,k,XWavelet,1, j,k,XWavelet,1)) << endl << endl;
                }
            }
        }
    }

    return 0;
}


template<Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Interval,Cons>& basis, int j, XType type)
{
    if (type==XBSpline) {
        return basis.mra.rangeI(j);
    }
    else {
        return basis.rangeJ(j);
    }
}

template<Construction Cons>
Range<int>
getRange(const Basis<T,Orthogonal,Interval,Cons>& basis, int j, XType type)
{
    if (type==XBSpline) {
        return basis.mra.rangeI(j);
    }
    else {
        return basis.rangeJ(j);
    }
}

template<Construction Cons>
Range<int>
getRange(const Basis<T,Orthogonal,R,Cons>& basis, int j, XType type)
{
    if (type==XBSpline) {
        return Range<int>(-15,15);
    }
    else {
        return Range<int>(-15,15);
    }
}

template<Construction Cons>
Range<int>
getRange(const Basis<T,Dual,Interval,Cons>& basis, int j, XType type)
{
    if (type==XBSpline) {
        //cout << basis.mra_.R_Left.firstCol() << " " << basis.mra_.R_Left.lastCol() << " : " << basis.mra_.R_Left << std::endl;
        //cout << basis.mra_.R_Right.firstCol() << " " << basis.mra_.R_Right.lastCol() << " : " << basis.mra_.R_Right << std::endl;
        //cout << basis.mra_.rangeI_L(j) << " " << basis.mra_.rangeI_I(j) << " " << basis.mra_.rangeI_R(j) << endl;
        return basis.mra_.rangeI_(j);
    }
    else {
        return basis.rangeJ_(j);
    }
}

