#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

/// FLENS typedefs: only required for the set up of right-hand side vectors
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;

///  Wavelet basis over an interval: Here we are exclusively (!!) considering $L_2$-orthonormal
///  multiwavelets.
typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;
typedef PrimalBasis::RefinementBasis                                  RefinementBasis;
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis>   Basis3D;

///  Definition of the (tensor product) wavelet preconditioner
typedef OptimizedH1Preconditioner3D<T,Basis3D>                        Preconditioner;

///  Underlying univariate bilinear form: As we are using $L_2$-orthonormal multiwavelets and
///  and are considering Poissons problem, we only need the univariate bilinear form
///  $a(v,w) = \int_0^1 v'(x) w'(x) dx$
typedef RefinementBasis::LaplaceOperator1D                            RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>        LaplaceOp1D;

///  Local operator in 1d for the above bilinear form
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D,LaplaceOp1D>            LocalOp1D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{A} \otimes \vec{\textrm{Id}} \otimes \vec{\textrm{Id}}$
typedef UniDirectionalLocalOperator<Index3D,XOne,LocalOp1D,
                                            NotXOne,Index2D>          UniDirectionalLocalOpXOne3D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{\textrm{Id}} \otimes \vec{A} \otimes \vec{\textrm{Id}}$
typedef UniDirectionalLocalOperator<Index3D,XTwo,LocalOp1D,
                                            NotXTwo,Index2D>          UniDirectionalLocalOpXTwo3D;

///  Set up of the two-dimensional operator for the evaluation of a matrix of the form
///  $\vec{\textrm{Id}} \otimes \vec{\textrm{Id}} \otimes \vec{A}$
typedef UniDirectionalLocalOperator<Index3D,XThree,LocalOp1D,
                                            NotXThree,Index2D>        UniDirectionalLocalOpXThree3D;

///  Aggregation of the two above two-dimensional operator in one class.
typedef CompoundLocalOperator<Index3D, UniDirectionalLocalOpXOne3D,
                                       UniDirectionalLocalOpXTwo3D,
                                       UniDirectionalLocalOpXThree3D> CompoundLocalOperator3D;

///  Righthandsides definitions: Here we are only considering right-hand side functions that can be
///  separated w.r.t.\ the coordinate direction, i.e., $f(x_1,x_2) = f_1(x_1) \otimes f_2(x)
///  or a sum of such functions.
typedef RHSWithPeaks1D<T,PrimalBasis>                                 Rhs1D;
typedef AdaptiveSeparableRhs<T,Index3D,Rhs1D,Rhs1D,Rhs1D>             AdaptiveSeparableRhsIntegral3D;
typedef CompoundRhs<T,Index3D,AdaptiveSeparableRhsIntegral3D,
                              AdaptiveSeparableRhsIntegral3D,
                              AdaptiveSeparableRhsIntegral3D>         CompoundRhsIntegral3D;

///  Some iterators we require for postprocessing.
typedef IndexSet<Index1D>::const_iterator                             const_set1d_it;
typedef IndexSet<Index3D>::const_iterator                             const_set3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator             coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator             coeff3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator       const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator       const_coeff3d_it;

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

T test(T x) {   return 100.; }

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

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index3D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_u3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u3,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f3);

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
    const char* treeType = "sparsetree";    // "gradedtree";
    bool sparsetree = false;
    if (strcmp(treeType,"sparsetree")==0) sparsetree = true;

    Timer time;
    T eps   = 1e-2;

    ///  Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis3D basis3d(basis,basis,basis);

    /// Initialization of the operators
    LaplaceOp1D                     laplaceOp1D(basis);
    LocalOp1D                       localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne3D     uniDirectionalOpXOne3D(localOp1D);
    UniDirectionalLocalOpXTwo3D     uniDirectionalOpXTwo3D(localOp1D);
    UniDirectionalLocalOpXThree3D   uniDirectionalOpXThree3D(localOp1D);
    CompoundLocalOperator3D         localOp3D(uniDirectionalOpXOne3D,uniDirectionalOpXTwo3D,
                                              uniDirectionalOpXThree3D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis3d,1.,1.,1.,0.);

    /// Initialization of the right-hand side
    DenseVectorT sing_pts_x, sing_pts_y, sing_pts_z;
    DenseMatrixT no_deltas, deltas_x, deltas_y, deltas_z;
    int order = 20;
    if (example==2) {  int order = 4+2*d; }
    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u2,sing_pts_y), fct_f2(f2,sing_pts_y);
    Function<T>                    fct_u3(u3,sing_pts_z), fct_f3(f3,sing_pts_z);
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

    ///  Initialization of $\mathbf{f}_\varepsilon$ required for the $\textbf{APPLY}$-based residual
    ///  for sufficiently small $\varepsilon$ (steered by wavelet levels in the implementation below)
    Coefficients<Lexicographical,T,Index3D> f_eps(SIZEHASHINDEX2D);
    cout << "Setting up reference right-hand side f_eps..." << endl;
    setUp_f_eps<T>(example, basis, Prec, f_eps,
                   rhs_u1_data, rhs_u1, rhs_u2_data, rhs_u2, rhs_u3_data, rhs_u3,
                   rhs_f1_data, rhs_f1, rhs_f2_data, rhs_f2, rhs_f3_data, rhs_f3);

    cout << "... finished." << endl;

    //T tol = 3e-07;
    cout << "Norm of f_eps: " << f_eps.norm() << endl;
    T tol = 1.;

    ///  Output files
    stringstream residual_error_filename;
    residual_error_filename << "error_multitree_mw_awgm_poisson3d_" << example << "_"
                            << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                            << residualType << "_" << treeType<< ".dat";
    ofstream residual_error_file(residual_error_filename.str().c_str());

    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index3D> u2(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index3D> tmp(SIZEHASHINDEX2D);
    //Coefficients<Lexicographical,T,Index3D> r_multitree(SIZEHASHINDEX2D);

    for (int iter=1; iter<=100; ++iter) {
        ///  To reduce storage requirement, we use the same vector for the $\textbf{APPLY}$-based
        ///  residual and the multitree based residual.
        Coefficients<Lexicographical,T,Index3D> r_approx(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index3D> r_eps(SIZEHASHINDEX2D);

        ///  Read Galerkin solution from file
        stringstream coefffilename;
        coefffilename << "coeff3d/coeff_multitree_mw_awgm_poisson3d_" << example << "_"
                      << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                      << residualType << "_" << treeType << "__" << iter << ".dat";
        readCoefficientsFromFile(u, coefffilename.str().c_str());

        if (u.size()==0) break;

        /*
        bool useSupportCenter=true;
        stringstream scatterplotfilename;
        scatterplotfilename << "scattercoeff_multitree_mw_awgm_poisson3d_" << example << "_"
                      << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                      << residualType << "_" << treeType << "__" << iter;
        plotScatterCoeff(u, basis3d, scatterplotfilename.str().c_str(), useSupportCenter);
        std::cerr << "Please hit enter." << std::endl;
        continue;
        */

        Index3D maxIndex, maxWaveletIndex;
        int *jmax = new int[1];
        int arrayLength = 1;
        getLevelInfo(u, maxIndex, maxWaveletIndex, jmax, arrayLength);
        int J = -100;
        for (int i=0; i<arrayLength; ++i) {
           if (J<jmax[i]) J=jmax[i];
        }

        cout << "Size of u: " << u.size() << endl;

        ///  Compute $\textbf{APPLY}$-based reference residual
        localOp3D.apply(u,r_eps,Prec,eps);
        r_eps -= f_eps;
        cout << "Reference computation of r_eps has finished." << endl;
        T exact_residual = r_eps.norm();
        cout << "Size of r_eps = " << r_eps.size() << endl;

        ///  Compute multitree-based residual
        time.start();
        extendMultiTree(basis3d, u, r_approx, residualType, sparsetree);
        localOp3D.eval(u,r_approx,Prec);
        for (coeff3d_it it=r_approx.begin(); it!=r_approx.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        cout << "Size of multitree residual = " << r_approx.size() << endl;
        T new_residual_time = time.elapsed();
        T new_residual_norm = r_approx.norm();
        int new_residual_length = r_approx.size();
        r_eps -= r_approx;
        T new_residual_diff = r_eps.norm();
        r_eps += r_approx;
        //cout << "diff = " << r << endl;


        T apply_residual_norm = 0., apply_residual_norm1 = 0., apply_residual_norm2 = 0.;
        int apply_residual_length = 0, apply_residual_length1 = 0, apply_residual_length2 = 0;
        T apply_residual_diff = 0., apply_residual_diff1 = 0., apply_residual_diff2 = 0.;
        T apply_residual_time = 0., apply_residual_time1 = 0., apply_residual_time2 = 0.;

        ///  Compute $\textbf{APPLY}$-based residual that has approximately the same accuracy as the
        ///  the multitree based residual
        while(1) {
            r_approx.clear();
            time.start();
            cout << "   Apply started..." << endl;
            localOp3D.apply(u,r_approx,Prec,tol/2.);
            cout << "   ... finished, output size = " << r_approx.size() << endl;
            //r -= f_eps;
            r_approx -= THRESH(f_eps,tol/2.,true,true);
            time.stop();
            apply_residual_norm2 = r_approx.norm();
            apply_residual_length2 = r_approx.size();
            apply_residual_time2 = time.elapsed();
            r_eps -= r_approx;
            apply_residual_diff2 = r_eps.norm();
            r_eps += r_approx;
            cerr << "   DEBUG: tol = " << tol
                 << ", apply_residual_diff = " << apply_residual_diff2
                 << ", new_residual_diff = " << new_residual_diff
                 << ", apply_residual_time = " << time.elapsed() << endl;;
            if (apply_residual_diff2<new_residual_diff)  {
                if (apply_residual_time2 < apply_residual_time1) {
                    apply_residual_norm = apply_residual_norm2;
                    apply_residual_diff = apply_residual_diff2;
                    apply_residual_length = apply_residual_length2;
                    apply_residual_time = apply_residual_time2;
                }
                else {
                    apply_residual_norm = apply_residual_norm1;
                    apply_residual_diff = apply_residual_diff1;
                    apply_residual_length = apply_residual_length1;
                    apply_residual_time = apply_residual_time1;
                }
                break;
            }
            else {
                apply_residual_norm1 = apply_residual_norm2;
                apply_residual_diff1 = apply_residual_diff2;
                apply_residual_length1 = apply_residual_length2;
                apply_residual_time1 = apply_residual_time2;
                tol *= 0.9;
            }
        }

        cout << u.size() << " " << r_eps.norm() << ":" << endl;
        cout << "new residual:   " << new_residual_norm << " " << new_residual_diff << " " << new_residual_length << " " << new_residual_time << endl;
        //cout << "new residual2:  " << new_residual_norm2 << " " << new_residual_diff2 << " " << new_residual_length2 << " " << new_residual_time2 << endl;
        cout << "apply residual: " << apply_residual_norm << " " << apply_residual_diff << " " << apply_residual_length << " " << apply_residual_time << endl << endl;

        residual_error_file << u.size() << " " << exact_residual << " "
                            << new_residual_norm << " " << new_residual_diff << " "
                            << new_residual_length << " " << new_residual_time << " "
                            << apply_residual_norm << " " << apply_residual_diff << " "
                            << apply_residual_length << " " << apply_residual_time <<  endl;

        eps = min(eps, 0.1*tol);

    }
    return 0;
}

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index3D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_u3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u3,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f3)
{
    int j0 = basis.j0;



    if (example==2) {
        for (int k=basis.mra.rangeIL(j0).firstIndex(); k<=basis.mra.rangeIL(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_u3_data[index1d] = rhs_u3(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            rhs_f3_data[index1d] = rhs_f3(index1d);
        }
        for (int k=basis.mra.rangeIR(j0).firstIndex(); k<=basis.mra.rangeIR(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_u3_data[index1d] = rhs_u3(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            rhs_f3_data[index1d] = rhs_f3(index1d);
        }

        for (int j=j0; j<=20; ++j) {
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_u3_data[index1d] = rhs_u3(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
                rhs_f3_data[index1d] = rhs_f3(index1d);
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_u3_data[index1d] = rhs_u3(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
                rhs_f3_data[index1d] = rhs_f3(index1d);
            }
        }

        for (const_coeff1d_it it_x=rhs_u1_data.begin(); it_x!=rhs_u1_data.end(); ++it_x) {
            for (const_coeff1d_it it_y=rhs_u2_data.begin(); it_y!=rhs_u2_data.end(); ++it_y) {
                for (const_coeff1d_it it_z=rhs_u3_data.begin(); it_z!=rhs_u3_data.end(); ++it_z) {
                    Index3D index((*it_x).first,(*it_y).first,(*it_z).first);
                    T val =  (   rhs_f1_data[(*it_x).first] * (*it_y).second  * (*it_z).second
                              + (*it_x).second * rhs_f2_data[(*it_y).first] * (*it_z).second
                              + (*it_x).second * (*it_y).second * rhs_f3_data[(*it_z).first]
                              ) * Prec(index);
                    if (fabs(val)>1e-16) f_eps[index] = val;
                }
            }
        }
    }

    std::cerr << "#Supp f_eps = " << f_eps.size() << std::endl;
}


/*
        extendMultiTree(basis3d, u, r_tmp, 1, sparsetree);
        localOp3D.eval(u,r_tmp,Prec,1);
        r_approx += r_tmp;

        r_tmp.clear();
        r_tmp = u;
        r_tmp.setToZero();
        extendMultiTree(basis3d, u, r_tmp, 2);
        localOp3D.eval(u,r_tmp,Prec,2);
        r_approx += r_tmp;

        r_tmp.clear();
        r_tmp = u;
        r_tmp.setToZero();
        extendMultiTree(basis3d, u, r_tmp, 3);
        localOp3D.eval(u,r_tmp,Prec,3);
        r_approx += r_tmp;
 */

/*
 *     DenseVectorT sing_pts;
    DenseMatrixT no_deltas;
    Function<T>                    fct_u1(test,sing_pts), fct_u2(u2,sing_pts), fct_u3(u3,sing_pts);
    RHSWithPeaks1D<T,PrimalBasis>  test_rhs_u1(basis, fct_u1, no_deltas, 8);
    RHSWithPeaks1D<T,PrimalBasis>  test_rhs_u2(basis, fct_u2, no_deltas, 8);
    RHSWithPeaks1D<T,PrimalBasis>  test_rhs_u3(basis, fct_u3, no_deltas, 8);

    T val2 = (  test_rhs_u1((*it_x).first) * test_rhs_u2((*it_y).first)
                              * test_rhs_u3((*it_z).first) ) * Prec(index);
                    cout << (*it_x).first << " " << (*it_y).first << " " << (*it_z).first << ": " << val << " " << val2 << endl;
 *
 */
