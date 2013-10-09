#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

typedef long double T;

///  Wavelet basis over the realline
typedef Basis<T, Orthogonal, R, Multi>                              PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Tensor product wavelet basis in two dimensions
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

///  Underlying bilinear form for $a(v_1 \otimes v_2,w_1 \otimes w_2) := \Big(\int_0^1 v_1'(x) w_1'(x) dx \Big) \Big(\int_0^1 v_2(x) w_2(x) dx \Big)
///  + \Big(\int_0^1 v_1(x) w_1(x) dx \Big)\Big(\int_0^1 v_2'(x) w_2'(x) dx \Big)$
///  Observe that the first definition of a corresponding operator is the standard definition and is
///  required for the computation of reference values. The second operator is solely for refinement
///  B-Spline bases as here, we only require evaluation of type
typedef LaplaceOperator1D<T,PrimalBasis>                            BilinearForm_x;
typedef RefinementBasis::LaplaceOperator1D                          RefinementBilinearForm_x;
typedef LaplaceOperator1D<T,PrimalBasis>                            BilinearForm_y;
typedef RefinementBasis::LaplaceOperator1D                          RefinementBilinearForm_y;

///  Local operator in 1d: These are required as building blocks for the two-dimensional operator.
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_x>                   LocalOp1D_x;
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_y>                   LocalOp1D_y;

///  Local operator in 2d for multiwavelet discretization: Allows for the evaluation of operators
///  $\vec{A} \otimes \vec{\textrm{Id}}$ and $\vec{\textrm{Id}} \otimes \vec{A}$ (here exemplarily
///  for two space dimensions. The extension to higher space dimensions is straightforward (see
///  the implementation.
typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D_x,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D_y,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

///  Iterators: Required for the calculation of reference solutions
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

///  A given coefficient vector is filled with random numbers.
void
getRandomCoefficientVector(Coefficients<Lexicographical,T,Index2D> &coeff);

///  Reference computation of the application of $\vec{\textrm{Id}} \otimes \vec{A}$
void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv);

///  Reference computation of the application of $\vec{A} \otimes \vec{\textrm{Id}}$
void
refComputationAIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AIv);

int main (int argc, char *argv[]) {
    cout.precision(6);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d j0_1 j0_2 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d     = atoi(argv[1]);
    int j0_x  = atoi(argv[2]);
    int j0_y  = atoi(argv[3]);

    int J     = atoi(argv[4]);

    int numOfIter=J;

    ///  Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis_x(d, j0_x);           // For L2_orthonormal MW bases
    PrimalBasis basis_y(d, j0_y);           // For L2_orthonormal MW bases

    RefinementBasis &refinementbasis_x = basis_x.refinementbasis;
    RefinementBasis &refinementbasis_y = basis_y.refinementbasis;
    Basis2D basis2d(basis_x,basis_y);

    ///  Operator initialization for univariate operators
    BilinearForm_x    Bil_x(basis_x);
    BilinearForm_y    Bil_y(basis_y);
    LocalOp1D_x localOperator_x(basis_x,basis_x,refinementbasis_x.LaplaceOp1D);
    LocalOp1D_y localOperator_y(basis_y,basis_y,refinementbasis_y.LaplaceOp1D);

    ///  Operator initialization for the two dimensional operator.
    UniDirectionalLocalOpXOne2D localop2d_x(localOperator_x);
    UniDirectionalLocalOpXTwo2D localop2d_y(localOperator_y);

    ///  We choose (arbitrarily) a two index and complete it to a multitree.
    Index2D index2d_1(Index1D(J,5,XWavelet),Index1D(J,8,XWavelet));
    Index2D index2d_2(Index1D(J-1,20,XWavelet),Index1D(J+2,17,XWavelet));
    Coefficients<Lexicographical,T,Index2D> v;
    completeMultiTree(basis2d, index2d_1, v, 0, true);
    completeMultiTree(basis2d, index2d_2, v, 0, true);
    getRandomCoefficientVector(v);
    cout << "v = " << v << endl;

    Index2D index2d_3(Index1D(J+1,12,XWavelet),Index1D(J+1,17,XWavelet));
    Coefficients<Lexicographical,T,Index2D> Av, IAv, AIv, IAv_ref, AIv_ref;
    completeMultiTree(basis2d, index2d_1, Av, 0, true);
    completeMultiTree(basis2d, index2d_2, Av, 0, true);
    completeMultiTree(basis2d, index2d_3, Av, 0, true);
    Av.setToZero();

    cout << "Av = " << Av << endl;

    IAv_ref = Av;
    AIv_ref = Av;
    IAv = Av;
    AIv = Av;

    ///   Computation of reference solutions
    refComputationAIv(Bil_x, v, AIv_ref);
    refComputationIAv(Bil_y, v, IAv_ref);

    ///  Call of the first local operator in two dimensions
    localop2d_x.eval(v, AIv);

    std::cerr << "Call of localop2d_x.eval finished." << std::endl;

    ///  Call of the second local operator in two dimensions
    localop2d_y.eval(v, IAv);

    std::cerr << "Call of localop2d_y.eval finished." << std::endl;

    Coefficients<Lexicographical,T,Index2D> IAv_error, AIv_error;
    IAv_error = IAv_ref - IAv;
    AIv_error = AIv_ref - AIv;
    cout << "IAv_error: " << IAv_error.norm(2.) << endl;
    cout << "AIv_error: " << AIv_error.norm(2.) << endl;

}

void
getRandomCoefficientVector(Coefficients<Lexicographical,T,Index2D> &coeff)
{
    for (coeff2d_it it=coeff.begin(); it!=coeff.end(); ++it) {
        (*it).second = T(rand()) / T(RAND_MAX);
    }
}

void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv)
{
    for (coeff2d_it row=IAv.begin(); row!=IAv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                val +=  Bil_y(row_y,col_y) * (*col).second;
            }
        }
        (*row).second = val;
    }
    return;
}

void
refComputationAIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AIv)
{
    for (coeff2d_it row=AIv.begin(); row!=AIv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                val +=  Bil_x(row_x,col_x) * (*col).second;
            }
        }
        (*row).second = val;
    }
    return;
}

