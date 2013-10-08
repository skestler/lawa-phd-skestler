#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef long double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:

///  Wavelet basis over the realline (so far, underlying algorithms are only implemented for $L_2$
///  orthonormal wavelet bases).
typedef Basis<T, Orthogonal, R, Multi>                              PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Underlying bilinear form for $a(v,w) := \int_{\mathbb{R}} w'(x) v'(x) dx$ and
///  $a(v,w) := \int_{\mathbb{R}} w(x) v(x) dx$ with optimized version only for refinement B-Splines
typedef LaplaceOperator1D<T,PrimalBasis>                            BilinearForm;
typedef RefinementBasis::LaplaceOperator1D                          RefinementBilinearForm;
typedef LaplaceOperator1D<T,RefinementBasis>                        RefinementBilinearFormTest;
//typedef IdentityOperator1D<T,PrimalBasis>                           BilinearForm;
//typedef RefinementBasis::IdentityOperator1D                         RefinementBilinearForm;
//typedef IdentityOperator1D<T,RefinementBasis>                       RefinementBilinearFormTest;

///  Initialization of a local operator in 1d - analogous to the interval case. The realline only
///  enters in the construction of the basis
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm>                     LocOp1D;

typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;

///  Bound for translation indices
int K = 20;

/// Construction of a random tree with (possibly) random values.
void
constructRandomTree(const PrimalBasis &basis, int J, bool withRandomValues,
                    TreeCoefficients1D<T> &LambdaTree, bool sparsetree, int seed);

/// Computation of a reference value $\mathbf{L}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
void
computeEvalLRef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Lv_tree);

/// Computation of a reference value $\mathbf{U}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
void
computeEvalURef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Uv_tree);

/// Computation of a reference value $\mathbf{A}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
void
computeEvalARef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Av_tree);

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d j0 J seed" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d    = atoi(argv[1]);
    int j0   = atoi(argv[2]);
    int J    = atoi(argv[3]);
    bool sparsetree = true;

    int seed = atoi(argv[4]);
    Timer time;

    /// Basis initialization
    PrimalBasis basis(d, j0);           // For L2_orthonormal MW bases
    RefinementBasis &refinementbasis = basis.refinementbasis;

    /// Operator initialization
    BilinearForm                    Bil(basis);
    RefinementBilinearFormTest      RefinementBilTest(refinementbasis);
    LocOp1D localOperator1D(basis,basis,refinementbasis.LaplaceOp1D);
    //LocOp1D localOperator1D(basis,basis,refinementbasis.IdentityOp1D);
    //LocOp1D localOperator1D(basis,basis,RefinementBilTest);


    int j=j0+2;
    for (int k1=-K; k1<=K; ++k1) {
        for (int k2=-K; k2<=K; ++k2) {
            if (abs(k1-k2)>9) continue;
            T val1 = RefinementBilTest(XBSpline, j, k1, XBSpline, j, k2);
            T val2 = refinementbasis.LaplaceOp1D(XBSpline, j, k1, XBSpline, j, k2);
            //T val2 = refinementbasis.IdentityOp1D(XBSpline, j, k1, XBSpline, j, k2);
            //if (fabs(val1-val2)>1e-13) {
                cout << "[" << j << ", (" << k1 << "," << k2 << ")]: " <<  val1
                << " " << val2 << " " << pow2i<T>(j) << endl;
            //}
        }
    }

    stringstream ct_filename;
    ct_filename << "comptime_locOp1d_" << d << "_" << j0 << "_" << J << ".dat";
    ofstream ct_file(ct_filename.str().c_str());


    for (int j=j0; j<=J; ++j) {
        T time_evalA = 0., time_evalU = 0., time_evalL = 0.;
        TreeCoefficients1D<T> v_tree(389,basis.j0), Av_tree(389,basis.j0),
                              Uv_tree(COEFFBYLEVELSIZE,basis.j0), Lv_tree(COEFFBYLEVELSIZE,basis.j0);
        TreeCoefficients1D<T> Av_ref_tree(COEFFBYLEVELSIZE,basis.j0),
                              Uv_ref_tree(COEFFBYLEVELSIZE,basis.j0),
                              Lv_ref_tree(COEFFBYLEVELSIZE,basis.j0);

        /// We construct random tree based coefficients where the output vector is initialized
        /// with zero values.
        constructRandomTree(basis, j, true, v_tree, sparsetree, seed);
        constructRandomTree(basis, j+1, false, Av_tree, sparsetree, seed+37);
        Av_ref_tree = Av_tree;
        Uv_tree     = Av_tree;
        Uv_ref_tree = Av_tree;
        Lv_tree     = Av_tree;
        Lv_ref_tree = Av_tree;

        cout << "#supp v_tree = " << v_tree.size() << endl;
        cout << "#supp Av_tree = " << Av_tree.size() << endl;

        cout << " ********* Testing evalA **********" << endl;
        int N = Av_tree.size() + v_tree.size();
        cout << "N = " << N << endl;
        time.start();

        /// Computation of a reference value for
        /// $\mathbf{A}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
        computeEvalARef(Bil, basis, v_tree, Av_ref_tree);

        time.stop();
        cout << "Elapsed time for EvalARef: " << time.elapsed() << endl;
        time.start();

        /// Computation of $\mathbf{A}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$
        /// by $\textbf{EvalS}$. Observe that the function call is the same as in the interval case.
        localOperator1D.eval(v_tree, Av_tree,"A");

        time.stop();
        time_evalA = time.elapsed();
        cout << "Elapsed time for evalA: " << time_evalA << endl;
        //cout << "v_tree = " << v_tree << endl;
        //cout << "Av_tree_ref = " << Av_ref_tree << endl;
        //cout << "Av_tree = " << Av_tree << endl;
        Av_ref_tree -= Av_tree;
        //cout << "diff = " << Av_ref_tree << endl;
        cout << "Error norm: " << Av_ref_tree.norm(2.) << endl;
        cout << " **********************************" << endl << endl;

        cout << " ********* Testing evalU **********" << endl;
        cout << "N = " << Uv_tree.size() + v_tree.size() << endl;
        time.start();

        /// Computation of a reference value for
        /// $\mathbf{U}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
        computeEvalURef(Bil, basis, v_tree, Uv_ref_tree);

        time.stop();
        time.start();

        /// Computation of $\mathbf{U}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$
        /// by $\textbf{EvalU}$. Observe that the function call is the same as in the interval case.
        localOperator1D.eval(v_tree, Uv_tree,"U");

        time.stop();
        time_evalU = time.elapsed();
        cout << "Elapsed time for evalU: " << time_evalU << endl;
        //cout << "v_tree = " << v_tree << endl;
        //cout << "Uv_tree_ref = " << Uv_ref_tree << endl;
        //cout << "Uv_tree = " << Uv_tree << endl;
        Uv_ref_tree -= Uv_tree;
        cout << "Error norm: " << Uv_ref_tree.norm(2.) << endl;
        cout << " **********************************" << endl << endl;

        cout << " ********* Testing evalL **********" << endl;
        cout << "N = " << Lv_tree.size() + v_tree.size() << endl;
        time.start();

        /// Computation of a reference value for
        /// $\mathbf{L}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$.
        computeEvalLRef(Bil, basis, v_tree, Lv_ref_tree);

        time.stop();
        time.start();

        /// Computation of $\mathbf{L}_{\widetilde{\Lambda}\times \Lambda} \mathbf{v}_{\Lambda}$
        /// by $\textbf{EvalL}$. Observe that the function call is the same as in the interval case.
        localOperator1D.eval(v_tree, Lv_tree,"L");

        time.stop();
        time_evalL = time.elapsed();
        cout << "Elapsed time for evalL: " << time_evalL << endl;
        //cout << "v_tree = " << v_tree << endl;
        //cout << "Uv_tree_ref = " << Uv_ref_tree << endl;
        //cout << "Uv_tree = " << Uv_tree << endl;
        Lv_ref_tree -= Lv_tree;
        cout << "Error norm: " << Lv_ref_tree.norm(2.) << endl;
        cout << " **********************************" << endl << endl;

        ct_file << N << " " << time_evalA << " " << time_evalU << " " << time_evalL << endl;
    }
    ct_file.close();

    return 0;
}

void
constructRandomTree(const PrimalBasis &basis, int J, bool withRandomValues,
                    TreeCoefficients1D<T> &LambdaTree, bool sparsetree, int seed)
{
    srand ( seed );
    T val = 0.;
    /*
    for (int k=-2*K; k<=2*K; ++k) {
        val = withRandomValues ? (rand() / (T)RAND_MAX) : 0.;
        LambdaTree[0].map.operator[](k) = val;
    }

    for (int j=basis.j0; j<=J; ++j) {
        //for (int k=-K; k<=K; ++k) {
        int k = 0;
            val = withRandomValues ? (rand() / (T)RAND_MAX) : 0.;
            LambdaTree.bylevel[j-basis.j0+1].map.operator[](k) = val;
        //}
    }
    */
    //cout << "LambdaTree = " << LambdaTree << endl;

    for (int k=-K; k<=K; ++k) {
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        LambdaTree[0].map.operator[](k) = val;
    }
    for (int j=basis.j0; j<=J; ++j) {
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        int random_k1 = rand() % (4*K+1) + 1;
        LambdaTree[j-basis.j0+1].map.operator[](random_k1) = val;
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        int random_k2 = rand() % (4*K+1) + 1;
        LambdaTree[j-basis.j0+1].map.operator[](random_k2) = val;
    }

    cout << "Before: LambdaTree = " << LambdaTree << endl;

    Coefficients<Lexicographical,T,Index1D> tmpTree, copyTmpTree;
    fromTreeCoefficientsToCoefficients(LambdaTree, copyTmpTree);
    for (coeff1d_it it=copyTmpTree.begin(); it!=copyTmpTree.end(); ++it) {
        completeMultiTree(basis, (*it).first, tmpTree, sparsetree);
    }

    for (coeff1d_it it=tmpTree.begin(); it!=tmpTree.end(); ++it) {
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        (*it).second = val;
    }

    fromCoefficientsToTreeCoefficients(tmpTree, LambdaTree);


    cout << "After: LambdaTree = " << LambdaTree << endl;

    /*
    if (sparsetree) cout << "Option "sparsetree=true" not available." << endl;
    for (int i=J-basis.j0+1; i>=2; --i) {
        CoefficientsByLevel<T> *currentlevel;
        currentlevel = &(LambdaTree.bylevel[i]);
        int j=basis.j0+i-1;
        for (const_coeffbylevel_it it=(*currentlevel).map.begin(); it!=(*currentlevel).map.end(); ++it) {
            long k = (*it).first;
            long k_first = std::max((int)k / 2 - 30, basis.rangeJ(j-1).firstIndex());
            long k_last  = std::min((int)k / 2 + 30, basis.rangeJ(j-1).lastIndex());
            for (int k1=k_first; k1<=k_last; ++k1) {
                if (overlap(basis.psi.support(j,k),basis.psi.support(j-1,k1))>0) {
                    val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
                    LambdaTree[i-1].map.operator[](k1) = val;
                }
            }
        }
    }
    */
}

void
computeEvalARef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Av_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Av;
    fromTreeCoefficientsToCoefficients(v_tree, v);
    fromTreeCoefficientsToCoefficients(Av_tree, Av);

    if (    (flens::IsSame<Basis<T,Orthogonal,R,Multi>, PrimalBasis>::value)
         && (flens::IsSame<IdentityOperator1D<T, Basis<T,Orthogonal,R,Multi> >, BilinearForm>::value) ) {
        cout << "Using orthogonality..." << endl;
        for (coeff1d_it row=Av.begin(); row!=Av.end(); ++row) {
            T val = 0.L;
            for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
                if (    (*row).first.xtype==(*col).first.xtype && (*row).first.j == (*col).first.j
                     && (*row).first.k == (*col).first.k)
                    val += (*col).second;
            }
            (*row).second = val;
        }
    }
    else {
        cout << "Not using orthogonality..." << endl;
        for (coeff1d_it row=Av.begin(); row!=Av.end(); ++row) {
            T val = 0.L;
            for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
                val +=  Bil((*row).first,(*col).first) * (*col).second;
            }
            (*row).second = val;
        }
    }
    fromCoefficientsToTreeCoefficients(Av, Av_tree);
}

void
computeEvalURef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Uv_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Uv;
    fromTreeCoefficientsToCoefficients(v_tree, v);
    fromTreeCoefficientsToCoefficients(Uv_tree, Uv);

    if (    (flens::IsSame<Basis<T,Orthogonal,R,Multi>, PrimalBasis>::value)
             && (flens::IsSame<IdentityOperator1D<T, Basis<T,Orthogonal,R,Multi> >, BilinearForm>::value) ) {
        cout << "Using orthogonality..." << endl;
        for (coeff1d_it row=Uv.begin(); row!=Uv.end(); ++row) {
            T val = 0.L;
            for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
                if (    (*row).first.xtype==(*col).first.xtype && (*row).first.j == (*col).first.j
                     && (*row).first.k == (*col).first.k)
                    val += (*col).second;
            }
            (*row).second = val;
        }
    }
    else {
        cout << "Not using orthogonality..." << endl;
        for (coeff1d_it row=Uv.begin(); row!=Uv.end(); ++row) {
            long double val = 0.L;
            for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
                if ( ((*row).first.xtype==XBSpline) ||
                     (((*row).first.xtype==XWavelet && (*col).first.xtype==XWavelet
                                              && (*row).first.j<=(*col).first.j)) ) {
                    val +=  Bil((*row).first,(*col).first) * (*col).second;
                }
            }
            (*row).second = val;
        }
    }
    fromCoefficientsToTreeCoefficients(Uv, Uv_tree);
}

void
computeEvalLRef(const BilinearForm &Bil, const PrimalBasis &basis,

                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Lv_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Lv;
    fromTreeCoefficientsToCoefficients(v_tree, v);
    fromTreeCoefficientsToCoefficients(Lv_tree, Lv);

    if (    (flens::IsSame<Basis<T,Orthogonal,R,Multi>, PrimalBasis>::value)
                 && (flens::IsSame<IdentityOperator1D<T, Basis<T,Orthogonal,R,Multi> >, BilinearForm>::value) ) {
        return;
    }
    else {
        for (coeff1d_it row=Lv.begin(); row!=Lv.end(); ++row) {
            long double val = 0.L;
            for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
                if ( !( (   (*row).first.xtype==XBSpline) ||
                        ( ( (*row).first.xtype==XWavelet && (*col).first.xtype==XWavelet
                                              && (*row).first.j<=(*col).first.j) ) ) ) {
                    val +=  Bil((*row).first,(*col).first) * (*col).second;
                }
            }
            (*row).second = val;
        }
    }
    fromCoefficientsToTreeCoefficients(Lv, Lv_tree);
}

