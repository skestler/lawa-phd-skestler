#include <iostream>
#include <lawa/lawa.h>

#include <applications/finance/initialconditions/initialconditions.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

/* ************************************ */
/* *** Typedefs for financial model *** */
/* ************************************ */

T strike = 1.;
T maturity = 1.;
T weight1 = 0.5, weight2 = 0.5;


const OptionTypenD optiontype = BasketPut;
OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketPutOption2D<T> > PayoffIntegral;
typedef RHS2D<T, PayoffIntegral, NoPreconditioner<T, Index2D>  >            PayoffIntegralRHS;

/*
const OptionTypenD optiontype = SumOfPuts;
OptionParameters2D<T,SumOfPuts> optionparameters(strike, strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedSumOfPutsOption2D<T> > PayoffIntegral;
*/

const ProcessType2D  processtype  = BlackScholes2D;
//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T r = 0.;
T sigma1 = 0.3;
T sigma2 = 0.2;

T rho = 0.;
T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T s1  = sigma1*sigma1, s2  = sigma2*sigma2;
T    critical_line_x1 = 0.6;
bool critical_above_x1 = true;

/*
T rho = 0.3;
T u11 = 0.95171801008793943164, u12 = 0.30697366218334239729, u21 = -0.30697366218334239729, u22 = 0.95171801008793943164;
T s1  = sqrt(13./2.*(199.+5.*sqrt(949)))/500., s2  = sqrt(13./2.*(199.-5.*sqrt(949)))/500.;
T    critical_line_x1 = 0.6;
bool critical_above_x1 = true;
*/

ProcessParameters2D<T,BlackScholes2D>   processparameters(r, sigma1, sigma2, rho, u11, u12, u21, u22);

/* ********************************************* */
/* *** Typedefs for numerical discretization *** */
/* ********************************************* */

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef RefinementBasis::IdentityOperator1D                         RefinementIdentityOp1D;


///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D>                      LocalOp1D;
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementIdentityOp1D>                     LocalIdentityOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;
typedef UniDirectionalLocalOperator<Index2D,XOne,LocalIdentityOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalIdentityOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalIdentityOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalIdentityOpXTwo2D;

typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;
typedef CompoundLocalOperator<Index2D,
                              UniDirectionalLocalIdentityOpXOne2D,
                              UniDirectionalLocalIdentityOpXTwo2D>  CompoundLocalIdentityOperator2D;

typedef ThetaTimeStepLocalOperator<Index2D, CompoundLocalOperator2D,
                                   CompoundLocalIdentityOperator2D> ThetaTimeStepLocalOperator2D;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef ThetaTimeStepSeparableRHS<T,Index2D,
                                  AdaptiveSeparableRhsIntegral2D,
                                  ThetaTimeStepLocalOperator2D>     ThetaTimeStepRhs2d;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;

typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      ThetaTimeStepRhs2d,Preconditioner>            ThetaTimeStepMultiTreeAWGM2D;

typedef LocalWeighting2D<T, Basis2D>                                LocalWeightingInitCond2D;

typedef MultiTreeAWGM<Index2D,Basis2D,
                      ThetaTimeStepLocalOperator2D,
                      PayoffIntegralRHS,
                      NoPreconditioner<T,Index2D> >                 ApproxL2AWGM2D;

typedef ThetaSchemeAWGM<Index2D, ThetaTimeStepMultiTreeAWGM2D>      ThetaSchemeMultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;


T f_t(T t)       {  return 0.; }
T f_x(T x)       {  return 0.; }
T f_y(T y)       {  return 0.; }

T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2);

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters);

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters,
                   const std::map<std::pair<T,T>,T> &refprices);

void
computeReferencePrice(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                      T inner_left1, T inner_right1, T inner_left2, T inner_right2, T h1, T h2,
                      const Coefficients<Lexicographical,T,Index2D> &u, int j,
                      Option2D<T,optiontype> &option2d,
                      ProcessParameters2D<T,BlackScholes2D> &processparameters);



void
readReferencePrice(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   T inner_left1, T inner_right1, T inner_left2, T inner_right2, T h1, T h2,
                   const Coefficients<Lexicographical,T,Index2D> &u, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters,
                   std::map<std::pair<T,T>,T> &refprices);

template <typename TruncatedOptionPayoff>
T
computeLinftyErrorInitCond(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                           const Coefficients<Lexicographical,T,Index2D> &u,T delta,
                           TruncatedOptionPayoff &truncatedpayoff2d);

void
plotComputedSolution(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                     const Coefficients<Lexicographical,T,Index2D> &u);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }

    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);
    T alpha = 0.4;
    T gamma = 0.025;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;
    int weightType = 1;
    size_t hashMapSize = 196613;

    T left_x1 = -2., right_x1 = 2.;
    T left_x2 = -2., right_x2 = 2.;
    T delta = 0.05;

    T theta = 0.5;
    T timestep_eps = 1e-2;
    int maxiterations =  50;  T init_cgtol = 1e-9;   // use maxiterations = 1 for "pure" sparse grid computation
    int numOfTimesteps = 128;
    T timestep = maturity/numOfTimesteps;

    int numOfMCRuns = 100000;

    int order = 4;

    bool useRefPrices = true;

    Timer time;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    DenseMatrixT U(2,2), tU(2,2), Q(2,2), QtU(2,2), UQtU(2,2);
    U  = u11, u12, u21, u22;
    tU = u11, u21, u12, u22;
    Q  = sigma1*sigma1, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2*sigma2;
    QtU  = Q(1,1)*tU(1,1)+Q(1,2)*tU(2,1), Q(1,1)*tU(1,2)+Q(1,2)*tU(2,2),
           Q(2,1)*tU(1,1)+Q(2,2)*tU(2,1), Q(2,1)*tU(1,2)+Q(2,2)*tU(2,2);
    UQtU = U(1,1)*QtU(1,1)+U(1,2)*QtU(2,1), U(1,1)*QtU(1,2)+U(1,2)*QtU(2,2),
           U(2,1)*QtU(1,1)+U(2,2)*QtU(2,1), U(2,1)*QtU(1,2)+U(2,2)*QtU(2,2);
    cout << "U Q U^T " << UQtU << endl;
    cout << "s1 = " << s1 << ", s2 = " << s2 << endl;

    T a1 = 0.5*s1/((right_x1-left_x1)*(right_x1-left_x1));
    T a2 = 0.5*s2/((right_x1-left_x1)*(right_x1-left_x1));

    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D, a1);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D, a2);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    LocalIdentityOp1D                    localIdentityOp1D(basis,basis,refinementbasis.IdentityOp1D);
    UniDirectionalLocalIdentityOpXOne2D  uniDirectionalIdentityOpXOne2D(localIdentityOp1D, 1.);
    UniDirectionalLocalIdentityOpXTwo2D  uniDirectionalIdentityOpXTwo2D(localIdentityOp1D, 0.);
    CompoundLocalIdentityOperator2D      localIdentityOp2D(uniDirectionalIdentityOpXOne2D,
                                                           uniDirectionalIdentityOpXTwo2D);
    ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,localOp2D,localIdentityOp2D);

    /// Initialization of preconditioner
    NoPreconditioner<T, Index2D> NoPrec;
    Preconditioner  Prec(basis2d, a1, a2, 1.);

    /// Initialization of integrals for initial condition and rhs
    DenseVectorT sing_pts_t, sing_pts_x(5), sing_pts_y(5);
    sing_pts_x = 0.1, 0.2, 0.3, 0.4, 0.5;
    sing_pts_y =  0.1, 0.2, 0.3, 0.4, 0.5;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    Function<T>                    fct_f_t(f_t,sing_pts_t);
    Function<T>                    fct_f_x(f_x,sing_pts_x), fct_f_y(f_y,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_x(basis, fct_f_x, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_y(basis, fct_f_y, no_deltas, order);

    Coefficients<Lexicographical,T,Index1D> rhs_f_x_data(SIZEHASHINDEX1D),
                                            rhs_f_y_data(SIZEHASHINDEX1D);

    AdaptiveSeparableRhsIntegral2D F_rhs(rhs_f_x, rhs_f_x_data, rhs_f_y, rhs_f_y_data);
    ThetaTimeStepRhs2d thetatimestep_F(fct_f_t,F_rhs,localThetaTimeStepOp2D);

    /// Initialization of integrals for initial condition and rhs
    Option2D<T,optiontype>         option2d(optionparameters);
    option2d.setNumberOfMCRuns(numOfMCRuns);

    TruncatedBasketPutOption2D<T> truncatedoption2d;
    //TruncatedSumOfPutsOption2D<T> truncatedoption2d;
    truncatedoption2d.setOption(option2d);
    truncatedoption2d.setTransformation(u11, u21, u12, u22);
    truncatedoption2d.setTruncation(left_x1, right_x1, left_x2, right_x2, 0, 0.2, 200.);
    truncatedoption2d.setCriticalLine_x1(critical_line_x1, critical_above_x1);

    PayoffIntegral payoffIntegral(basis2d, truncatedoption2d,
                                  left_x1, right_x1, left_x2, right_x2, true, 0.05, order);
    PayoffIntegralRHS payoffIntegralRHS(payoffIntegral, NoPrec);

    Coefficients<Lexicographical,T,Index2D> u(hashMapSize), u0(hashMapSize), f(hashMapSize);

    std::map<std::pair<T,T>,T> refprices;
    if (useRefPrices) {
        readReferencePrice(basis2d, left_x1, right_x1, left_x2, right_x2,
                           -0.1, 0.1, -0.1, 0.1, 0.02, 0.02, u, 0, option2d, processparameters,
                           refprices);
    }

    std::stringstream filename;

    if (optiontype == BasketPut) {
        if (useRefPrices) {
            filename << "basketputoption2d_conv2_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                     << "_lx2_" << left_x2 << "_rx2_" << right_x2
                     << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                     << processparameters << ".txt";
        }
        else {
            filename << "basketputoption2d_conv_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                     << "_lx2_" << left_x2 << "_rx2_" << right_x2
                     << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                     << processparameters << ".txt";
        }
    }
    else if (optiontype == SumOfPuts) {
        if (useRefPrices) {
            filename << "sumofputsoption2d_conv2_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                     << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
                     << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                     << processparameters << ".txt";
        }
        else {
            filename << "sumofputsoption2d_conv_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                     << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
                     << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                     << processparameters << ".txt";
        }
    }
    else {
        cerr << "Option type does not exist." << endl;
        exit(1);
    }
    std::ofstream convfile(filename.str().c_str());

    for (int j=0; j<=J; ++j) {
        T timestep_eps = 1e-10;
        int maxL2Iterations = 18; u.clear();
        //int maxL2Iterations = 12+j; u.clear();
        //int maxL2Iterations = 15; u.clear();
        getSparseGridVector(basis2d, u, j, (T)0.);
    //for (timestep_eps=0.1; timestep_eps>=1e-6; timestep_eps*=0.5) {
        //int maxL2Iterations = 100; u.clear();
        //getSparseGridVector(basis2d, u, 0, (T)0.);

        LocalWeightingInitCond2D localWeightingInitCond2D;
        localWeightingInitCond2D.setDomain(left_x1,right_x1,left_x2,right_x2);
        localWeightingInitCond2D.setBasis(&basis2d);
        localWeightingInitCond2D.setWeightType(weightType);

        ApproxL2AWGM2D approxL2_solver(basis2d, localThetaTimeStepOp2D, payoffIntegralRHS, NoPrec);
        approxL2_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);
        approxL2_solver.approxL2(u, timestep_eps, localWeightingInitCond2D.weight, maxL2Iterations);
        cout << "Approximation of initial condition finished." << endl;

        T maxerror_initcond = computeLinftyErrorInitCond(basis2d, left_x1, right_x1,
                                                         left_x2, right_x2,
                                                          u, 0.25, truncatedoption2d);
        cout << "Maxerror initcond = " << maxerror_initcond << endl;

        /// Initialization of multi tree based adaptive wavelet Galerkin method
        ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
                                                              thetatimestep_F, Prec);
        thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false,
                                           hashMapSize);


        int strategy = 2;
        ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver);
        thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
                                  init_cgtol, strategy);

        int avDof = 0, maxDof = 0., terminalDof;
        thetascheme.solve(u, avDof, maxDof, terminalDof, j);
        cerr << "Computation of u has finished." << endl;

        plotComputedSolution(basis2d, left_x1, right_x1, left_x2, right_x2, u);

        T maxerror = 0., maxerror1 = 0., maxerror2 = 0.;
        if (useRefPrices) {
            maxerror1 = computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2,
                                          u,delta,-1,option2d, processparameters, refprices);
            maxerror2 = computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2,
                                          u,delta,-1,option2d, processparameters);
            convfile << timestep << " " << -1 << " " << avDof << " " << maxDof << " " << terminalDof << " "
                     << maxerror1 << " " << maxerror2 << " " << maxerror_initcond << " "
                     << timestep_eps << " " << numOfMCRuns << " " << delta << endl;
        }
        else {
            maxerror = computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2,
                                          u,delta,-1,option2d, processparameters);
            convfile << timestep << " " << -1 << " " << avDof << " "
                     << maxerror << " " << maxerror_initcond << " "
                     << timestep_eps << " " << numOfMCRuns << " " << delta << endl;
        }
        cerr << "Computation of errors has finished." << endl;
    }
    return 0;

}

T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2)
{
    T RightmLeft_x1 = right_x1-left_x1, SqrtRightmLeft_x1 = std::sqrt(right_x1-left_x1);
    T RightmLeft_x2 = right_x2-left_x2, SqrtRightmLeft_x2 = std::sqrt(right_x2-left_x2);

    T ret = 0.;
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        int   j1 = (*it).first.index1.j,     j2 = (*it).first.index2.j;
        int   k1 = (*it).first.index1.k,     k2 = (*it).first.index2.k;
        XType e1 = (*it).first.index1.xtype, e2 = (*it).first.index2.xtype;

        T val_x1 = (1./SqrtRightmLeft_x1) * basis2d.first.generator(e1).operator()((x1-left_x1)/(RightmLeft_x1),j1,k1,0);
        T val_x2 = (1./SqrtRightmLeft_x2) * basis2d.second.generator(e2).operator()((x2-left_x2)/(RightmLeft_x2),j2,k2,0);

        ret += (*it).second * val_x1 * val_x2;
    }
    return ret;
}

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters)
{
    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "bs2d_basketput_awgm_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "bs2d_sumofputs_awgm_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

    T maxerror = 0.;
    T h1 = (delta*right_x1-delta*left_x1)/10.;
    T h2 = (delta*right_x2-delta*left_x2)/10.;
    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=delta*left_x1; x1<=delta*right_x1; x1+=h1) {
        //for (T x2=left_x2; x2<=right_x2; x2+=0.03125) {
        for (T x2=delta*left_x2; x2<=delta*right_x2; x2+=h2) {
            /*
            T S1 = strike*std::exp(u11*x1+u21*x2+(0.5*sigma1*sigma1-r)*maturity);
            T S2 = strike*std::exp(u12*x1+u22*x2+(0.5*sigma2*sigma2-r)*maturity);
            T exact = std::exp(r*maturity)*option2d.value(processparameters,S1,S2,0);
            T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
            T approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1, x2);
            */
            T S1    = strike*std::exp(x1);
            T S2    = strike*std::exp(x2);
            T exact = option2d.value(processparameters,S1,S2,0);
            T x1hat = u11*x1+u12*x2-u11*(0.5*sigma1*sigma1-r)*maturity-u12*(0.5*sigma2*sigma2-r)*maturity;
            T x2hat = u21*x1+u22*x2-u21*(0.5*sigma1*sigma1-r)*maturity-u22*(0.5*sigma2*sigma2-r)*maturity;
            T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
            T approx =std::exp(-r*maturity)*evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat, x2hat);
            plotfile << x1 << " " << x2 << " " << exact << " " << approx << " " << payoff << endl;
            maxerror = std::max(maxerror, fabs(approx-exact));
        }
        plotfile << endl;
    }
    plotfile.close();

    return maxerror;
}

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters,
                   const std::map<std::pair<T,T>,T> &refprices)
{
    typedef std::map<std::pair<T,T>,T>::const_iterator const_map_it;

    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "bs2d_basketput2_awgm_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "bs2d_sumofputs2_awgm_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

    T maxerror = 0.;
    for (const_map_it it=refprices.begin(); it!=refprices.end(); ++it) {
        T x1    = (*it).first.first;
        T x2    = (*it).first.second;
        T exact = (*it).second;

        //T approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1, x2);
        T x1hat = u11*x1+u12*x2-u11*(0.5*sigma1*sigma1-r)*maturity-u12*(0.5*sigma2*sigma2-r)*maturity;
        T x2hat = u21*x1+u22*x2-u21*(0.5*sigma1*sigma1-r)*maturity-u22*(0.5*sigma2*sigma2-r)*maturity;
        T approx =std::exp(-r*maturity)*evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat, x2hat);

        T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
        maxerror = std::max(maxerror, fabs(approx-exact));
        plotfile << x1 << " " << x2 << " " << exact << " " << approx << " " << payoff << endl;
    }
    plotfile.close();

    return maxerror;
}


void
computeReferencePrice(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                      T inner_left1, T inner_right1, T inner_left2, T inner_right2, T h1, T h2,
                      const Coefficients<Lexicographical,T,Index2D> &u, int j,
                      Option2D<T,optiontype> &option2d,
                      ProcessParameters2D<T,BlackScholes2D> &processparameters)
{
    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "bs2d_refprices_basketput_awgm_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "bs2d_refprices_sumofputs_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

    for (T x1=inner_left1; x1<=inner_right1; x1+=h1) {
        for (T x2=inner_left2; x2<=inner_right2; x2+=h2) {
            T S1 = strike*std::exp(x1+(0.5*sigma1*sigma1-r)*maturity);
            T S2 = strike*std::exp(x2+(0.5*sigma2*sigma2-r)*maturity);
            T exact = std::exp(r*maturity)*option2d.value(processparameters,S1,S2,0);
            T approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1, x2);
            plotfile << x1 << " " << x2 << " " << exact << " " << approx << endl;
        }
        plotfile << endl;
    }
    plotfile.close();
}

void
readReferencePrice(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   T inner_left1, T inner_right1, T inner_left2, T inner_right2, T h1, T h2,
                   const Coefficients<Lexicographical,T,Index2D> &u, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,BlackScholes2D> &processparameters,
                   std::map<std::pair<T,T>,T> &refprices)
{
    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "refprices/bs2d_refprices_basketput_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "refprices/bs2d_refprices_sumofputs_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::cerr << "Try to open file: " << filename.str().c_str() << std::endl;
    std::ifstream infile (filename.str().c_str());
    if (infile.is_open()) {
        cout << "File is open, ready to read..." << endl;
    }
    else {
        cout << "File is not open. Exit..." << endl; exit(1);
    }

    std::string line;
    while(std::getline( infile, line, '\n' )) {
        std::string field1, field2, field3, field4;
        std::istringstream line_ss(line);
        std::getline( line_ss, field1, ' ' );
        std::getline( line_ss, field2, ' ' );
        std::getline( line_ss, field3, ' ' );
        std::getline( line_ss, field4, ' ' );

        T x1     = atof(field1.c_str());
        T x2     = atof(field2.c_str());
        T exact  = atof(field3.c_str());
        T approx = atof(field4.c_str());
        cerr.precision(16);
        if (exact!=0) {
            refprices[std::pair<T,T>(x1,x2)] = approx;
            cerr.precision(16);
            cerr << "Using approx values: " << x1 << " " << x2 << " " << approx << endl;
        }
    }
}

template <typename TruncatedOptionPayoff>
T
computeLinftyErrorInitCond(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                           const Coefficients<Lexicographical,T,Index2D> &u,T delta,
                           TruncatedOptionPayoff &truncatedpayoff2d)
{
    std::stringstream coeffsfilename;
    coeffsfilename << "coeffs_initcond_" << u.size();
    plotScatterCoeff(u, basis2d, coeffsfilename.str().c_str(), true, left_x1, right_x1, left_x2, right_x2);

    std::stringstream filename;
    filename << "debug_initcond_" << u.size() << ".txt";
    ofstream file(filename.str().c_str());
    T maxerror = 0.;
    //T h1 = (delta*right_x1-delta*left_x1)/200.;
    //T h2 = (delta*right_x2-delta*left_x2)/200.;
    T h1 = (right_x1-left_x1)/100.;
    T h2 = (right_x2-left_x2)/100.;
    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=left_x1; x1<=right_x1; x1+=h1) {
        //for (T x2=left_x2; x2<=right_x2; x2+=0.03125) {
        for (T x2=left_x2; x2<=right_x2; x2+=h2) {
            T exact = truncatedpayoff2d.payoff(x1,x2);
            T approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1, x2);
            file << x1 << " " << x2 << " " << exact << " " << approx << endl;

            if (x1>=delta*left_x1 && x1<=delta*right_x1 && x2>=delta*left_x2 && x2<=delta*right_x2) {
                maxerror = std::max(maxerror, fabs(exact-approx));
            }

        }
        file << endl;
    }
    file.close();
    return maxerror;
}

void
plotComputedSolution(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                     const Coefficients<Lexicographical,T,Index2D> &u)
{
    std::stringstream coeffsfilename;
    coeffsfilename << "coeffs_endpoint_" << u.size();
    plotScatterCoeff(u, basis2d, coeffsfilename.str().c_str(), true, left_x1, right_x1, left_x2, right_x2);

    std::stringstream filename;
    filename << "solution_" << u.size() << ".txt";
    ofstream file(filename.str().c_str());
    T maxerror = 0.;
    T h1 = (right_x1-left_x1)/100.;
    T h2 = (right_x2-left_x2)/100.;
    for (T x1=left_x1; x1<=right_x1; x1+=h1) {
        for (T x2=left_x2; x2<=right_x2; x2+=h2) {
            T approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1, x2);
            file << x1 << " " << x2 << " " << approx << endl;
        }
        file << endl;
    }
    file.close();
}

