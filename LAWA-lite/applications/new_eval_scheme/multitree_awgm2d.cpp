#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;


typedef double T;

/// FLENS typedefs: only required for the set up of right-hand side vectors
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

///  Wavelet basis over an interval
typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
//typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;

///  Tensor product wavelet basis in two dimensions
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;


typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>    WeightedPDEBilinearForm;
//typedef AdaptiveWeightedPDEOperator1D<T,Orthogonal,Interval,MultiRefinement>    WeightedPDEBilinearForm;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying local operators the evaluation of $a(v,w) := \Big(\int_0^1 p(x) v(x) w(x) dx \Big)$
///  or $a(v,w) := \Big(\int_0^1 p(x) v'(x) w'(x) dx \Big)$
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        WeightedPDEBilinearForm>                    LocalWeightedPDEOp1D;

///  Local operator in 2d: Build from the above univariate components
typedef LocalOperator2D<LocalWeightedPDEOp1D,
                        LocalWeightedPDEOp1D>                       LocalWeightedPDEOp2D;

///  Aggregation of two or more local operators in one class.
typedef CompoundLocalOperator<Index2D, LocalWeightedPDEOp2D,
                              LocalWeightedPDEOp2D>                 CompoundLocalOperator2D;

///  Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D >                                  SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;

typedef RHS<T,Index2D,SumOfSeparableRhsIntegral2D,
            Preconditioner>                                         SumOfSeparableRhs;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

///  Routine that writes a (multitree) index set to a file (for debugging)
void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, const char *name, int example, int d, T threshTol, int ell, int nr);

///  Routine that reads a (multitree) index set from a file (for debugging)
void
readIndexSetFromFile(IndexSet<Index2D> &Lambda,  int example, int d, T threshTol, int ell, int nr);

///  Routine for the set up of an initial sparse grid index set
void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

///  Routine for the extension of the right-hand side index set (for testing).
void
extendRHSIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j);

///  Construct multitree extension of an existing multitree (see Section 6.5)
void
simpleExtendMultiTree(const Basis2D &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda);

///  Wrapper for the matrix vector multiplication permitting the (easy) measurement of computation
///  times.
void
mv(CompoundLocalOperator2D &localOperator2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time);

T p1(T x)  {   return (x-0.5)*(x-0.5)+1.; /*1.;*/  }

T dp1(T x) {   return 2*(x-0.5);          /*0.;*/  }

T p2(T y)  {   return (y-0.5)*(y-0.5)+1.; /*1.;*/  }

T dp2(T y) {   return 2*(y-0.5);          /*0.;*/  }

/*
int example = 1;

T u1(T x)
{
    if      (0<=x && x<1./3.)    return 10*exp(x)-10.;
    else if (1./3<x && x<=2./3.) return -10.+10.*exp(1./3.)+(-10.+(30.+90.*(-(2./3.)+x))*(-(1./3.)+x))*(-(1./3.)+x);
    else                         return 10*exp(-(x-1.))-10.;
}

T u2(T y)
{
    if      (0<=y && y<1./3.)    return 10*exp(y)-10.;
    else if (1./3<y && y<=2./3.) return -10.+10.*exp(1./3.)+(-10.+(30.+90.*(-(2./3.)+y))*(-(1./3.)+y))*(-(1./3.)+y);
    else                         return 10*exp(-(y-1.))-10.;
}

T du1(T x)
{
    if      (0<=x && x<1./3.)    return 10*exp(x);
    else if (1./3<x && x<=2./3.) return 20-180.*x+270.*x*x;
    else                         return -10*exp(-(x-1.));
}

T du2(T y)
{
    if      (0<=y && y<1./3.)    return 10*exp(y);
    else if (1./3<y && y<=2./3.) return 20-180.*y+270.*y*y;
    else                         return -10*exp(-(y-1.));
}

T ddu1(T x)
{
    if      (0<=x && x<1./3.)   return 10*exp(x);
    else if (1./3<x && x<2./3.) return -180. + 540*x;
    else                        return 10*exp(-(x-1.));
}

T ddu2(T y)
{
    if      (0<=y && y<1./3.)   return 10*exp(y);
    else if (1./3<y && y<2./3.) return -180. + 540*y;
    else                        return 10*exp(-(y-1.));
}
T L2norm_x_sq = 10./567.*(23011.-26775.*exp(1./3.)+7560.*exp(2./3.));
T L2norm_y_sq = L2norm_x_sq;
T H1semi_x_sq = 20./3.*(-11.+15.*exp(2./3.));
T H1semi_y_sq = H1semi_x_sq;
T H1seminorm_squared = L2norm_x_sq*H1semi_y_sq + H1semi_x_sq*L2norm_y_sq;

T p_u1(T x) {   return p1(x)*u1(x); }

T f1(T x)   {   return -p1(x)*ddu1(x)-dp1(x)*du1(x); }

T p_u2(T y) {   return p2(y)*u2(y); }

T f2(T y)   {   return -p2(y)*ddu2(y)-dp2(y)*du2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }
*/

int example = 2;

T u1(T x)   {    return 1.; }

T u2(T y)   {    return 1.; }

T du1(T x)  {    return 0.; }

T du2(T y)  {    return 0.; }

T ddu1(T x) {    return -10.;   }

T ddu2(T y) {    return -10.;   }

//T H1seminorm_squared = 14.05770149526849;
T H1seminorm_squared = 11.68737041271871;

T p_u1(T x) {   return u1(x); }

T f1(T x)   {   return -ddu1(x); }

T p_u2(T y) {   return u2(y); }

T f2(T y)   {   return -ddu2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }

/*
int example = 3;

T u1(T x)   {    return x*x*(1-x)*(1-x); }

T u2(T y)   {    return y*y*(1-y)*(1-y); }

T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }

T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }

T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }

T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }


T H1seminorm_squared = 2.*(13./630. * 1./616.);

T p_u1(T x) {   return p1(x)*u1(x); }

T f1(T x)   {   return -p1(x)*ddu1(x)-dp1(x)*du1(x); }

T p_u2(T y) {   return p2(y)*u2(y); }

T f2(T y)   {   return -p2(y)*ddu2(y)-dp2(y)*du2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }
*/


int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);
    bool adaptive=true;
    T threshTol = 0.6;
    T r_norm = 0.1;
    T gamma = 0.2;
    int ell=1;
    Timer time;

    /// Basis initialization
    PrimalBasis       basis(d,d_,j0);
    //PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis,basis);

    /// Initialization of operator with weighted coefficients
    DenseVectorT p1_singPts, p2_singPts;
    Function<T> reaction_coeff(p1, p1_singPts);
    Function<T> convection_coeff(p1, p1_singPts);
    Function<T> diffusion_coeff(p1, p1_singPts);
    WeightedPDEBilinearForm       WeightedLaplaceBil( basis.refinementbasis,reaction_coeff,convection_coeff,diffusion_coeff,10,true,true,false);
    WeightedPDEBilinearForm       WeightedIdentityBil(basis.refinementbasis,reaction_coeff,convection_coeff,diffusion_coeff,10,false,true,true);

    /// Initialization of local operator
    LocalWeightedPDEOp1D          localLaplaceOp1D( basis, basis, WeightedLaplaceBil);
    LocalWeightedPDEOp1D          localIdentityOp1D(basis, basis, WeightedIdentityBil);
    LocalWeightedPDEOp2D          localLaplaceIdentityOp2D(localLaplaceOp1D,localIdentityOp1D);
    LocalWeightedPDEOp2D          localIdentityLaplaceOp2D(localIdentityOp1D,localLaplaceOp1D);
    localLaplaceIdentityOp2D.setJ(7);
    localIdentityLaplaceOp2D.setJ(7);
    CompoundLocalOperator2D       localOperator2D(localLaplaceIdentityOp2D,localIdentityLaplaceOp2D);

    /// Initialization of preconditioner
    //HelmholtzBilinearForm2D  HelmholtzBil2D(basis2d,0.);
    //Preconditioner           Prec(HelmholtzBil2D);
    Preconditioner           Prec(basis2d,1.,1.,0.);


    DenseVectorT sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    if (example==1) {
        sing_pts_x.engine().resize(2); sing_pts_x(1) = 1./3.; sing_pts_x(2) = 2./3.;
        sing_pts_y.engine().resize(2); sing_pts_y(1) = 1./3.; sing_pts_y(2) = 2./3.;
        deltas_x.engine().resize(2,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 10.+10.*exp(1./3.);
                                       deltas_x(2,1) = 2./3.; deltas_x(2,2) = 20.+10.*exp(1./3.);
        deltas_y.engine().resize(2,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = 10.+10.*exp(1./3.);
                                       deltas_y(2,1) = 2./3.; deltas_y(2,2) = 20.+10.*exp(1./3.);
    }
    SeparableFunction2D<T> SepFunc1(f1, sing_pts_x, p_u2, sing_pts_y);
    SeparableFunction2D<T> SepFunc2(p_u1, sing_pts_x, f2, sing_pts_y);
    int order = 40;
    if (example==2) order=3*d;

    ///  Initialization of the right-hand side exploiting the product structure
    SeparableRhsIntegral2D      rhsintegral_x(basis2d, SepFunc1, deltas_x, no_deltas, order);
    SeparableRhsIntegral2D      rhsintegral_y(basis2d, SepFunc2, no_deltas, deltas_y, order);
    SumOfSeparableRhsIntegral2D rhsintegral2d(rhsintegral_x,rhsintegral_y);
    SumOfSeparableRhs           F(rhsintegral2d,Prec);

    ///  Initialization of hash map vectors: observe that the size of `SIZEHASHINDEX2D` has a
    ///  strong effect on the computation times. These quantities are defined in
    ///  `lawa/methods/adaptive/datastructures/index.h`. They determine the number of initial buckets
    ///  in the hash map. If you choose them too small, rehashing will lead to non-optimal computation
    ///  times. If you choose them too large, the algorithm is slowned down. As a rule of thumb, you
    ///  may choose them 3-4 times larger then the number of degrees of freedom you want to use.
    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> f(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> r(SIZELARGEHASHINDEX2D),
                                            p(SIZEHASHINDEX2D),
                                            Ap(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> P(SIZEHASHINDEX2D);

    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis,Lambda,0,gamma);
    //readIndexSetFromFile(Lambda,example,d,threshTol,1,13);

    stringstream filename;
    filename << "multitree_awgm2d_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << ".txt";
    ofstream file(filename.str().c_str());
    file.precision(16);

    for (int iter=0; iter<=100; ++iter) {

        //readIndexSetFromFile(Lambda,example,d,threshTol,1,iter);
        //writeIndexSetToFile(Lambda,"Lambda",example,d,threshTol,ell,iter);

        if (Lambda.size()>400000) break;
        cout << endl;

        cout << "******** Iteration " << iter << endl;

        cout << "   Current size of Lambda " << Lambda.size() << endl;
        f = F(Lambda);

        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(Lambda,Lambda_x,Lambda_y);
        int jmin_x=0, jmax_x=0, jmin_y=0, jmax_y=0;
        getMinAndMaxLevel(Lambda_x,jmin_x,jmax_x);
        getMinAndMaxLevel(Lambda_y,jmin_y,jmax_y);
        cout << "   Max level = (" << jmax_x << ", " << jmax_y << ")" << endl;

        r.clear();
        p.clear();
        Ap.clear();
        FillWithZeros(Lambda,r);
        FillWithZeros(Lambda,p);
        FillWithZeros(Lambda,Ap);

        int maxIterations = 200;
        T tol;
        if (adaptive) tol=std::min((T)1e-2,1e-2*r_norm);
        else          tol=1e-8;

        T alpha, beta, rNormSquare, rNormSquarePrev;
        T dummy=0.;
        cerr << "   Computing preconditioner." << endl;

        ///  We precompute and store preconditioner values for better performance.
        for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
            if (P.find((*it))==P.end()) P[(*it)] = Prec(*it);
        }
        cerr << "   Computing matrix vector product for initial residual." << endl;

        ///  Implementation of the cg algorithm
        mv(localOperator2D, P, Lambda, u, r, dummy);
        r -= f;
        p -= r;
        rNormSquare = r*r;
        int cg_iters=0;
        T mv_time=0.;
        for (cg_iters=0; cg_iters<maxIterations; ++cg_iters) {
            if (sqrt(rNormSquare)<=tol) {
                cerr << "      CG stopped with error " << sqrt(rNormSquare) << endl;
                break;
            }
            T time=0.;
            mv(localOperator2D, P, Lambda, p, Ap, time);
            mv_time += time;
            T pAp = p * Ap;
            alpha = rNormSquare/pAp;
            u += alpha*p;
            r += alpha*Ap;
            rNormSquarePrev = rNormSquare;
            rNormSquare = r*r;
            cerr << "      Current error in cg: " << std::sqrt(rNormSquare) << endl;
            beta = rNormSquare/rNormSquarePrev;
            p *= beta;
            p -= r;
        }

        mv(localOperator2D, P, Lambda, u, Ap, dummy);
        T uAu = Ap*u;
        T fu  = f*u;
        T H1error =  sqrt(fabs(H1seminorm_squared - 2*fu + uAu));
        cerr.precision(20);
        std::cerr << "---> H1-error: " << H1error << ": " << H1seminorm_squared << " " << fu  << " " << uAu << endl;
        cerr.precision(6);

        /*
        stringstream plotfilename;
        plotfilename << "multitree_awgm2d_Linfty_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << iter;
        std::cerr << "Start plot with u.size() == " << u.size() << std::endl;
        plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0.3, 0.7, 0.3, 0.7, 0.1, plotfilename.str().c_str());
        std::cerr << "Finished plot with u.size() == " << u.size() << std::endl;
        */
        /*
        Coefficients<Lexicographical,T,Index1D> tmp;
        for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
            if ((*it).first.index1.xtype==XBSpline && (*it).first.index1.k==2) tmp[(*it).first.index2] = (*it).second;
        }
        cout << "tmp = " << tmp << endl;
        stringstream coefffilename;
        coefffilename << "u_" << iter;
        plotCoeff(tmp, basis, coefffilename.str().c_str(), false, true);
        */

        ///  Computation of the multitree extension of the current Galerkin index set: Here, we
        ///  provide two different possibilities where the second one adds further indices at the
        ///  boundary.
        IndexSet<Index2D> checkLambda = Lambda;
        Timer time;
        if (example!=2) {
            IndexSet<Index2D> C_Lambda = Lambda;
            for (int i=1; i<=ell; ++i) {
                C_Lambda = C(C_Lambda, (T)1., basis2d, true);
                for (const_set2d_it it=C_Lambda.begin(); it!=C_Lambda.end(); ++it) {
                    simpleExtendMultiTree(basis2d,(*it),checkLambda);
                }
            }
        }
        else {
            IndexSet<Index2D> LambdaBoundary;
            extendRHSIndexSet(basis, LambdaBoundary, std::max(jmax_x,jmax_y)+1);
            //f = F(LambdaBoundary);
            IndexSet<Index2D> C_Lambda = Lambda;
            C_Lambda += LambdaBoundary;
            for (int i=1; i<=ell; ++i) {
                C_Lambda = C(C_Lambda, (T)1., basis2d, true);
                time.start();
                for (const_set2d_it it=C_Lambda.begin(); it!=C_Lambda.end(); ++it) {
                    simpleExtendMultiTree(basis2d,(*it),checkLambda);
                }
                time.stop();
            }
        }
        time.start();
        for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
            if (P.find((*it))==P.end()) P[(*it)] = Prec(*it);
        }
        time.stop();
        cout << "   Elapsed time for preconditioner: " << time.elapsed() << endl;

        //writeIndexSetToFile(checkLambda,"checkLambda",example,d,threshTol,ell,iter);
        cerr << "   Extension of rhs finished." << endl;
        T time_res;

        ///   Computation of the approximate residual vector
        f = F(checkLambda);
        mv(localOperator2D, P, checkLambda, u, r, time_res);
        time_res=0.;
        mv(localOperator2D, P, checkLambda, u, r, time_res); // only for checking computation times!
        r-=f;

        cerr << "   Matrix vector for residual computation finished." << endl;

        r_norm = r.norm(2.);
        file << Lambda.size() << " " << checkLambda.size()
             << " " << H1error << " " << r_norm << " " << mv_time/cg_iters << " " << time_res << endl;
        cerr << "---> H1-error = " << H1error << ", res = " << r_norm << endl;

        if (adaptive) {
            T P_Lambda_r_norm_square = 0.;
            for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
                P_Lambda_r_norm_square += std::pow(r[(*it)],(T)2.);
                r.erase((*it));
            }

            ///  Bulk chasing with subsequent multitree completion
            r = THRESH(r,threshTol*r.norm(2.));
            for (const_coeff2d_it it=r.begin(); it!=r.end(); ++it) {
                simpleExtendMultiTree(basis2d,(*it).first,Lambda);
            }

        }
        else {
            getSparseGridIndexSet(basis,Lambda,iter+1,gamma);
        }
        cerr << "   Computation of new Lambda finished." << endl;
        r.clear();
    }

    return 0;
}

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, const char *name, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << name << "_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << nr << ".dat";
    ofstream file(filename.str().c_str());
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        file << *it << endl;
    }
    file.close();
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "indexsets/Lambda2d_" << example << "_" << d << "_"
             << threshTol << "_" << ell << "_" << nr << ".dat";
    std::ifstream infile (filename.str().c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
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
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0+i2-1;
            for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
                Lambda.insert(Index2D(col,row));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0+i1-1;
        for (int i2=1; i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) {
                continue;
            }
            int j2=j0+i2-1;
            for (long k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

void
extendRHSIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int J)
{
    std::cerr << "   extendRHSIndexSet adds level up to " << J << std::endl;
    IndexSet<Index1D> LambdaBoundary;
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaBoundary.insert(Index1D(basis.j0,k,XBSpline));
    }
    for (int j=basis.j0; j<=J; ++j) {
        for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
            LambdaBoundary.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
            LambdaBoundary.insert(Index1D(j,k,XWavelet));
        }
    }
    for (const_set1d_it it_x=LambdaBoundary.begin(); it_x!=LambdaBoundary.end(); ++it_x) {
        for (const_set1d_it it_y=LambdaBoundary.begin(); it_y!=LambdaBoundary.end(); ++it_y) {
            Lambda.insert(Index2D(*it_x,*it_y));
        }
    }
}

void
simpleExtendMultiTree(const Basis2D &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int offset = 7;

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {   return;                 }
    else                                    {   Lambda.insert(index2d); }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    Support<T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (index_x.j==j0_x) {
        for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
            Support<T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
                if (Lambda.find(new_index2d)==Lambda.end()) simpleExtendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.first.rangeJ(index_x.j-1).firstIndex(),index_x.k/2 - offset);
        long k_last  = std::min((long)basis.first.rangeJ(index_x.j-1).lastIndex(),index_x.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
                if (Lambda.find(new_index2d)==Lambda.end()) simpleExtendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }

    Support<T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (index_y.j==j0_y) {
        for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
            Support<T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
                if (Lambda.find(new_index2d)==Lambda.end()) simpleExtendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.second.rangeJ(index_y.j-1).firstIndex(),index_y.k/2 - offset);
        long k_last  = std::min((long)basis.second.rangeJ(index_y.j-1).lastIndex(),index_y.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
                if (Lambda.find(new_index2d)==Lambda.end()) simpleExtendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    return;
}

void
mv(CompoundLocalOperator2D &localOperator2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time)
{
    Timer timer;
    Av.setToZero();
    FillWithZeros(Lambda,Av);


    cout << "      MV start." << endl;
    timer.start();
    localOperator2D.eval(v,Av,P);
    timer.stop();
    cout << "      MV stop." << endl;
    time = timer.elapsed();


    std::cerr << "      MV: dof = " << Av.size() << ", time = " << time  << std::endl;
}

/*
T P_Lambda_r_norm_square = 0.;
for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
    P_Lambda_r_norm_square += std::pow(r[(*it)],2.);
    r.erase((*it));
}
T bulkalpha = 0.6;
T threshbound = std::sqrt(1-bulkalpha*bulkalpha) * r.norm(2.)/std::sqrt(T(r.size()));
Coefficients<Bucket,T,Index2D> r_bucket;
std::cerr << "      norm of r = " << r_norm << std::endl;
std::cerr << "      size of r = " << r.size() << std::endl;

r_bucket.bucketsort(r, threshbound);
IndexSet<Index2D> refinements;
for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
    P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
    int addDoF = r_bucket.addBucketToIndexSet(refinements,i);
    std::cerr << "      Added " << addDoF << " indices, now ||P_{Lambda}r ||_2 = "
              << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << bulkalpha*r_norm << std::endl;
    if (P_Lambda_r_norm_square >= bulkalpha*r_norm*bulkalpha*r_norm) {
        int addDoF = r_bucket.addBucketToIndexSet(refinements,i+1);
        std::cerr << "      Added " << addDoF << " indices, now ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
        break;
    }
}

cout << " Refinements: " << refinements << endl;
for (const_set2d_it it=refinements.begin(); it!=refinements.end(); ++it) {
    extendMultiTree2(basis2d,(*it),offset,Lambda);
}
*/

/*
stringstream coeff_filename;
if (adaptive) {
    coeff_filename << "u_adap_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << u.size();
}
else {
    coeff_filename << "u_sg_" << example << "_" << d << "_" << gamma << "_" << ell << "_" << u.size();
}
plotScatterCoeff2D<T,Index2D,PrimalBasis,PrimalBasis>(u, basis, basis, coeff_filename.str().c_str());
*/

/*
Coefficients<Lexicographical,T,Index1D> u_multi;
for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
    XType xtype1 = (*it).index1.xtype;
    int j1 = (*it).index1.j;
    long k1 = (*it).index1.k;
    if (xtype1==XBSpline && k1==basis.mra.rangeI(basis.j0).firstIndex()) {
        u_multi[(*it).index2] = 1.;
    }
}
stringstream filename2;
filename2 << "u_multi_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << iter << ".dat";
plotCoeff<T,PrimalBasis>(u_multi, basis, filename2.str().c_str(), false, true);
*/

