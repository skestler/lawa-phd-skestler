#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema> Basis1D;
//typedef Basis<T,Orthogonal,Interval,Multi> Basis1D;


int main (int argc, char *argv[]) {

    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(1);
    }

    int d = atoi(argv[1]);      // polynomial order
    int d_ =atoi(argv[2]);      // vanishing moments
    int j0 = atoi(argv[3]);     // minimal level
    int J = atoi(argv[4]);      // maximum level

    Basis1D basis(d,d_,j0);
    //Basis1D basis(d,j0);
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();

    T a = 0., b = 1.;

    ofstream plotfile_scaling("scaling.txt");
    for (T x=a; x<=b; x+=pow2i<T>(-6)) {
        plotfile_scaling << x;
        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            T val = basis.generator(XBSpline)(x,j0,k,0); // = basis.mra.phi(x,j,k,0);
            plotfile_scaling << " " << val;
        }
        plotfile_scaling << endl;
    }

    ofstream plotfile_wavelet("wavelet.txt");
    for (T x=a; x<=b; x+=pow2i<T>(-6)) {
        plotfile_wavelet << x;
        for (int j=j0; j<=J-1; ++j) {
            for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
                T val = basis.generator(XWavelet)(x,j,k,0);  // = basis.psi(x,j,k,0);
                plotfile_wavelet << " " << val;
            }
        }
        plotfile_wavelet << endl;
    }

    return 0;
}
