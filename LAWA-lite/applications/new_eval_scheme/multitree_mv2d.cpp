/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:

///  Wavelet basis over an interval
typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

///  Underlying bilinear form
typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>    BilinearForm_x;
typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>    RefinementBilinearForm_x;
typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>    BilinearForm_y;
typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>    RefinementBilinearForm_y;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_x>                   LocalOp1D_x;
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_y>                   LocalOp1D_y;
typedef LocalOperator2D<LocalOp1D_x, LocalOp1D_y>                   LocalOp2D;

///  Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, const char* indexset, int example, int d,
                     T threshTol, int ell, int nr);

void
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
                           Coefficients<Lexicographical,T,Index2D> &coeff);

void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv);

void
refComputationLIIAv(BilinearForm_x &Bil_y, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv);

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv);
void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv);

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv);

T p1(T x)  {   return (x-0.5)*(x-0.5)+1.; /*1.;*/  }

T dp1(T x) {   return 2*(x-0.5);          /*0.;*/  }

T p2(T y)  {   return (y-0.5)*(y-0.5)+1.; /*1.;*/  }

T dp2(T y) {   return 2*(y-0.5);         /*0.;*/  }


int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif
    cout.precision(6);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    int numOfIter=J;
    bool calcRefSol=false;

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    DenseVectorT p1_singPts, p2_singPts;
    Function<T> reaction_coeff(p1, p1_singPts);
    Function<T> convection_coeff(p1, p1_singPts);
    Function<T> diffusion_coeff(p1, p1_singPts);
    BilinearForm_x  RefinementBil_x(basis.refinementbasis,reaction_coeff,convection_coeff,diffusion_coeff,10,true,true,false);
    BilinearForm_y  RefinementBil_y(basis.refinementbasis,reaction_coeff,convection_coeff,diffusion_coeff,10,false,true,true);
    BilinearForm_x  Bil_x(basis,reaction_coeff,convection_coeff,diffusion_coeff,10,true,true,false);
    BilinearForm_y  Bil_y(basis,reaction_coeff,convection_coeff,diffusion_coeff,10,false,true,true);
    LocalOp1D_x localOperator_x(basis,basis,RefinementBil_x);
    LocalOp1D_y localOperator_y(basis,basis,RefinementBil_y);

    LocalOp2D   localop2d(localOperator_x,localOperator_y);
    localop2d.setJ(9);

    Timer time;

    ofstream file2("multitree_mv2d.dat");

    T old_time = 1.;
    T old_N = 1.;
    T time_evalAA1 = 0.;
    T time_intermediate1=0., time_intermediate2=0.,
                  time_IAv1=0., time_IAv2=0., time_LIv=0., time_UIv=0.;
    T time_intermediate1_old=0., time_intermediate2_old=0.,
      time_IAv1_old=0., time_IAv2_old=0., time_LIv_old=0., time_UIv_old=0.;
    int N = 0, N_old = 0;

    for (int j=0; j<=numOfIter; ++j) {

        IndexSet<Index2D> checkLambda, Lambda;


        T threshTol = 0.6;
        int ell=1;
        int example = 2;
        readIndexSetFromFile(Lambda,"Lambda",example,d,threshTol,ell,j);
        readIndexSetFromFile(checkLambda,"checkLambda",example,d,threshTol,ell,j);
        cout << "Size of Lambda:      " << Lambda.size() << endl;
        cout << "Size of checkLambda: " << checkLambda.size() << endl;
        if (Lambda.size()==0) return 0;

        Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D);


        getRandomCoefficientVector(Lambda,v);

        if (calcRefSol) {
            T time_evalAA1 = 0.;
            Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
            Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);
            Coefficients<Lexicographical,T,Index2D> IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref;
            IndexSet<Index1D> checkLambda_x;
            for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
                checkLambda_x.insert((*it).index1);
                LIIAv[*it] = 0.;
                LIIAv_ref[*it] = 0.;
                IAUIv[*it] = 0.;
                IAUIv_ref[*it] = 0.;
                AAv_ref[*it] = 0.;
            }
            for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
                Index1D col_x = (*col).first.index1;
                Index1D col_y = (*col).first.index2;
                for (const_set2d_it row=checkLambda.begin(); row!=checkLambda.end(); ++row) {
                    Index1D row_x = (*row).index1;
                    Index1D row_y = (*row).index2;
                    if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                          || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                        Support<T> col_supp_x = basis.generator(col_x.xtype).support(col_x.j,col_x.k);
                        Support<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                        if (overlap(col_supp_x,row_supp_x)>0) {
                            Index2D index(col_x,row_y);
                            IAv_ref[index] = 0.;
                        }
                    }
                }
            }
            for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
                Index1D col_x = (*col).first.index1;
                Index1D col_y = (*col).first.index2;
                for (const_set1d_it row=checkLambda_x.begin(); row!=checkLambda_x.end(); ++row) {
                    Index1D row_x = (*row);
                    if (     (row_x.xtype==XBSpline)
                          || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j <= col_x.j)) {
                        Support<T> col_supp_x = basis.generator(col_x.xtype).support(col_x.j,col_x.k);
                        Support<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                        if (overlap(col_supp_x,row_supp_x)>0) {
                            Index2D index(row_x,col_y);
                            UIv_ref[index] = 0.;
                        }
                    }
                }
            }
            cout << "Size of checkLambda: " << checkLambda.size() << endl;
            cout << "Size of Lambda:      " << Lambda.size() << endl;
            cout << "Size of IAv:         " << IAv_ref.size() << endl;
            cout << "Size of UIv:         " << UIv_ref.size() << endl;
            cout << "Size of v:           " << v.size() << endl;

            cout << "Reference calculation started..." << endl;
            refComputationIAv(Bil_y, v, IAv_ref);
            cout << "IAv_ref finished." << endl;
            refComputationLIIAv(Bil_x, IAv_ref, LIIAv_ref);
            cout << "LIIAv_ref finished." << endl;
            refComputationUIv(Bil_x, v, UIv_ref);
            cout << "UIv_ref finished." << endl;
            refComputationIAUIv(Bil_y, UIv_ref, IAUIv_ref);
            cout << "IAUIv_ref finished." << endl;
            refComputationAAv(Bil_x,Bil_y, v, AAv_ref);
            cout << "AAv_ref finished." << endl;
            cout << "Reference calculation finished." << endl;
            cout << "New scheme started..." << endl;
            time.start();
            localop2d.debug_eval(v, LIIAv, IAUIv, IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref);
            time.stop();
            time_evalAA1 = time.elapsed();
            cout << "New scheme finished." << endl;
        }
        else {
            Coefficients<Lexicographical,T,Index2D> AAv(SIZEHASHINDEX2D);

            for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
                AAv[*it] = 0.;
            }
            N = v.size() + checkLambda.size();
            cout << "**** New scheme started ****" << endl;
            cout << "   #v = " << Lambda.size() << endl;

            localop2d.eval(v, AAv, time_intermediate1, time_intermediate2,
                           time_IAv1, time_IAv2, time_LIv, time_UIv);

            AAv.setToZero();
            time.start();
            localop2d.eval(v, AAv, time_intermediate1, time_intermediate2,
                           time_IAv1, time_IAv2, time_LIv, time_UIv);
            time.stop();
            time_evalAA1 = time.elapsed();
            cout << "   N = " << N << ", time = " << time_evalAA1 << " -> ratio new / old = "
                 << (T)v.size()/old_N << ", " << time_evalAA1/old_time
                 << ", msec/dof = " << 1000.*time_evalAA1/N << endl;
            cout << "   " << N << " " << time_intermediate1 << " " <<  time_intermediate2 << " " << time_IAv1
                 << " " << time_IAv2 << " " << time_LIv << " " << time_UIv << endl;
            cout << "   " << T(N)/N_old << " : " << time_intermediate1/time_intermediate1_old
                               << " " << time_intermediate2/time_intermediate2_old
                               << " " << time_IAv1/time_IAv1_old << " " << time_IAv2/time_IAv2_old
                               << " " << time_LIv/time_LIv_old << " " << time_UIv/time_UIv_old << endl;
            cout << "**** New scheme finished ****" << endl << endl;
            // Attention: For large output sets, computation times are not exactly linear since
            // hash maps are too small.
            N_old = N;
            time_intermediate1_old=time_intermediate1; time_intermediate2_old=time_intermediate2;
            time_IAv1_old=time_IAv1; time_IAv2_old=time_IAv2; time_LIv_old=time_LIv;
            time_UIv_old=time_UIv;
        }
        file2 << v.size() << " " << checkLambda.size() << " " << time_evalAA1 << endl;
        old_N = v.size();
        old_time = time_evalAA1;
    }


    return 0;
}



void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, const char* indexset, int example, int d,
                     T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "indexsets/" << indexset << "_" << example << "_" << d << "_"
             << threshTol << "_" << ell << "_" << nr << ".dat";
    std::ifstream infile (filename.str().c_str());
    if (!infile.is_open()) {
        cerr << "   Indexset file " << filename.str().c_str()  << " is not open." << endl;
    }

    std::string line;
    std::string field1, field2, field3, field4, field5, field6;
    while(std::getline( infile, line, '\n' )) {
        std::istringstream line_ss(line);
        std::getline( line_ss, field1, ',' );
        std::getline( line_ss, field2, ',' );
        std::getline( line_ss, field3, ',' );
        std::getline( line_ss, field4, ',' );
        std::getline( line_ss, field5, ',' );
        std::getline( line_ss, field6, ',' );
        int j1,j2;
        long k1,k2;

        j1 = atoi(field2.c_str());
        k1 = atol(field3.c_str());
        j2 = atoi(field5.c_str());
        k2 = atol(field6.c_str());

        if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Got " << field1 << ", could not read file." << std::endl;
            exit(1); return;
        }
    }
}

void
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
                           Coefficients<Lexicographical,T,Index2D> &coeff)
{
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        coeff[*it] = T(rand()) / T(RAND_MAX);
    }
    return;
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
refComputationLIIAv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv)
{
    for (coeff2d_it row=LIIAv.begin(); row!=LIIAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=IAv.begin(); col!=IAv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                      || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                    val +=   Bil_x(row_x,col_x) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
}

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv)
{
    for (coeff2d_it row=UIv.begin(); row!=UIv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (    (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet
                                 && row_x.j<=col_x.j)) ) {
                    val += Bil_x(row_x,col_x) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
    return;
}

void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv)
{
    for (coeff2d_it row=IAUIv.begin(); row!=IAUIv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=UIv.begin(); col!=UIv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                val +=   Bil_y(row_y,col_y) * (*col).second;
            }
        }
        (*row).second = val;
    }
}

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv)
{
    for (coeff2d_it row=AAv.begin(); row!=AAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
                val +=   Bil_x(row_x,col_x) * Bil_y(row_y,col_y) * (*col).second;
        }
        (*row).second = val;
    }
}
