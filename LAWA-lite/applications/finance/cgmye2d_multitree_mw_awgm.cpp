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

typedef Basis<T,Orthogonal,R,Multi>                                 PrimalBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

/* ************************************ */
/* *** Typedefs for financial model *** */
/* ************************************ */


T strike = 1.;
T maturity = 1.;
T weight1 = 0.5, weight2 = 0.5;
/*
const OptionTypenD optiontype = BasketPut;
OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketPutOption2D<T> > PayoffIntegral;
typedef RHS2D<T, PayoffIntegral, NoPreconditioner<T, Index2D>  >            PayoffIntegralRHS;
*/

const OptionTypenD optiontype = SumOfPuts;
OptionParameters2D<T,SumOfPuts> optionparameters(strike, strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedSumOfPutsOption2D<T> > PayoffIntegral;
typedef RHS2D<T, PayoffIntegral, NoPreconditioner<T, Index2D>  >            PayoffIntegralRHS;

const ProcessType2D  processtype  = CGMYeUnivariateJump2D;
//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T r = 0.;
T sigma1 = 0.3;
T sigma2 = 0.2;
T rho = 0.;
//T k_C1 = 1., k_G1 = 7.4, k_M1 = 8.5, k_Y1 = 0.8;
//T k_C2 = 1., k_G2 = 6.5, k_M2 = 9.5, k_Y2 = 1.1;
T k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 1.25;
T k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;

//T k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 0.2;
//T k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;
//T k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 1.55;
//T k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;

T    critical_line_x1 = 0.6;
bool critical_above_x1 = true;

ProcessParameters2D<T,CGMYeUnivariateJump2D>   processparameters(r, sigma1, sigma2, rho,
                                                                 k_C1,  k_G1, k_M1, k_Y1,
                                                                 k_C2,  k_G2, k_M2, k_Y2);


/* ********************************************* */
/* *** Typedefs for numerical discretization *** */
/* ********************************************* */

//typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;
typedef FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>           CGMYeOp2D;
typedef ThetaTimeStepLocalOperator<Index2D,CGMYeOp2D>               ThetaTimeStepLocalOperator2D;
typedef DiagonalMatrixPreconditioner2D<T,Basis2D,CGMYeOp2D>         Preconditioner;

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

template <typename TruncatedOption>
void
PlotInitialCondition(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                     const Coefficients<Lexicographical,T,Index2D> &u,
                     TruncatedOption &truncatedOption, T &maxError, T &innerMaxError);

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int N,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters);

void
getSparseGridIndexSet(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                      Coefficients<Lexicographical,T,Index2D> &u, int j);

void
getLeftAndRightTranslationIndices(const PrimalBasis &basis, T left_x, T right_x,
                                  int j, XType type, long &left_k, long &right_k);

int main (int argc, char *argv[]) {

    cout.precision(16);
    if (argc!=8) {
        cout << "Usage: " << argv[0] << " d j0 J R1_1 R2_1 R1_2 R2_2" << endl;
        return 0;
    }

    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J   = atoi(argv[3]);
    T alpha = 0.2;
    T gamma = 0.025;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;
    int weightType = 1;
    size_t hashMapSize = 196613;
    T R1_1 = atof(argv[4]);
    T R2_1 = atof(argv[5]);
    T left_x1 = -R1_1, right_x1 = R2_1;
    T R1_2 = atof(argv[6]);
    T R2_2 = atof(argv[7]);
    T left_x2 = -R1_2, right_x2 = R2_2;
    T delta = 0.05;

    T theta = 0.5;
    T timestep_eps = 1e-6;
    int maxiterations =  1;  T init_cgtol = 1e-9;   // use maxiterations = 1 for "pure" sparse grid computation
    int numOfTimesteps = 16;
    T timestep = maturity/numOfTimesteps;
    int maxL2Iterations = 1;

    int numOfMCRuns = 100000;

    int order = 4;

    bool useRefPrices = false;

    Timer time;

    /// Basis initialization
    PrimalBasis     basis(d,j0);
    Basis2D         basis2d(basis,basis);

    cout << "Process parameters: " << processparameters << endl;
    CGMYeOp2D                    cgmyeOp2D(basis2d, processparameters,
                                           0., 1., 0., 1., 10);
    ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,cgmyeOp2D);

    /// Initialization of preconditioner
    NoPreconditioner<T, Index2D> NoPrec;
    //Preconditioner  Prec(basis2d, sigma1*sigma1, sigma2*sigma2, 1.);
    Preconditioner  Prec(cgmyeOp2D);


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

    ///   Initialization of the truncated payoff function (see, e.g., p. 183 for a visualization)
    TruncatedSumOfPutsOption2D<T> truncatedoption2d;
    //TruncatedBasketPutOption2D<T> truncatedoption2d;
    truncatedoption2d.setOption(option2d);
    truncatedoption2d.setTruncation(left_x1, right_x1, left_x2, right_x2, 0, 0.1, 100.);
    truncatedoption2d.setCriticalLine_x1(critical_line_x1, critical_above_x1);

    ///  Initialization of the integral for approximating the truncated payoff function
    PayoffIntegral payoffIntegral(basis2d, truncatedoption2d,
                                  0., 1., 0., 1., true, 0.2, order);
    PayoffIntegralRHS payoffIntegralRHS(payoffIntegral, NoPrec);

    Coefficients<Lexicographical,T,Index2D> u(hashMapSize), f(hashMapSize);

    std::stringstream filename;

    if (optiontype == BasketPut) {
        filename << "basketputoption2d_conv2_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                 << "_lx2_" << left_x2 << "_rx2_" << right_x2
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "sumofputsoption2d_conv2_awgm_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
                 << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    std::ofstream convfile(filename.str().c_str());


    for (int j=0; j<=J; ++j) {
        u.clear();
        std::cerr << "Computation of initial index set." << std::endl;
        getSparseGridIndexSet(basis2d, left_x1, right_x1, left_x2, right_x2, u, j);
        std::cerr << "Computation of adopted sparse grid index set finished. Size: " << u.size() << std::endl;

        ///  Initialization of the local weighting of the residual defined in Eq. (8.126). Here, a
        ///  value of zero refers to no weighting. However, for experimental purposes also other
        ///  values may be used.
        LocalWeightingInitCond2D localWeightingInitCond2D;
        localWeightingInitCond2D.setDomain(left_x1,right_x1,left_x2,right_x2);
        localWeightingInitCond2D.setBasis(&basis2d);
        localWeightingInitCond2D.setWeightType(0);

        ///  Setting the compression level for the two-dimensional CGMY operator (see Eq. (8.132))
        cgmyeOp2D.setCompressionLevel(j, j);

        ///  Initializing and calling the AWGM solver for approximating the initial condition
        ApproxL2AWGM2D approxL2_solver(basis2d, localThetaTimeStepOp2D, payoffIntegralRHS, NoPrec);
        approxL2_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);
        approxL2_solver.approxL2(u, timestep_eps, localWeightingInitCond2D.weight, maxL2Iterations);
        cout << "Approximation of initial condition finished." << endl;

        T maxError = 0., innerMaxError = 0.;
        //PlotInitialCondition(basis2d, left_x1, right_x1, left_x2, right_x2, u, truncatedoption2d, maxError, innerMaxError);
        cout << "Size of u: " << u.size() << " : " << maxError << " " << innerMaxError << endl;

        ///  Initializing the AWGM solver for each time-step
        ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
                                                          thetatimestep_F, Prec);
        thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false,
                                           hashMapSize);

        ///  For an adaptive strategy, i.e., an AWGM solve in each time-step (or, as an optimization only
        ///  only in certain time-steps), stratey should be set to 1. For using a fixed index set
        ///  that has a sparse grid structure but that remains fixed throughout the $\theta$-scheme,
        ///  strategy should be set to 0. For more details, please see the implementation of
        ///  the $\theta$-scheme in lawa/application/adaptive/solver/thetaschemeawgm.tcc.
        int strategy = 0;

        ///  Initializing and calling the $\theta$-scheme
        ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver);
        thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
                                  init_cgtol, strategy);

        int avDof = 0, maxDof = 0., terminalDof;
        thetascheme.solve(u, avDof, maxDof, terminalDof, 4);
        cerr << "Computation of u has finished." << endl;

        T maxErrorPrice = computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2, u,
                                             delta, avDof, option2d, processparameters);

        convfile << u.size() << " " << maxErrorPrice  << std::endl;
        cerr << "---> Max. error: " << maxErrorPrice << endl;
    }

    return 0;
}

T
evaluate(const Basis2D &basis2d, const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2)
{
    T ret = 0.;
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        int   j1 = (*it).first.index1.j,     j2 = (*it).first.index2.j;
        int   k1 = (*it).first.index1.k,     k2 = (*it).first.index2.k;
        XType e1 = (*it).first.index1.xtype, e2 = (*it).first.index2.xtype;

        T val_x1 = basis2d.first.generator(e1).operator()(x1,j1,k1,0);
        T val_x2 = basis2d.second.generator(e2).operator()(x2,j2,k2,0);

        ret += (*it).second * val_x1 * val_x2;
    }
    return ret;
}

template <typename TruncatedOption>
void
PlotInitialCondition(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                     const Coefficients<Lexicographical,T,Index2D> &u,
                     TruncatedOption &truncatedOption, T &maxError, T &innerMaxError)
{
    std::stringstream coeffsfilename;
    coeffsfilename << "cgmye2d_coeffs_initcond_" << u.size();
    plotScatterCoeff2D(u, basis2d.first, basis2d.second, coeffsfilename.str().c_str());


    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "cgmye2d_basketput_awgm_initcond_" << "_strike_" << strike
                 << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2 << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "cgmye2d_sumofputs_awgm_initcond_" << "_strike_" << strike
                 << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2 << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

    maxError = 0.;
    innerMaxError = 0.;
    T h1 = (right_x1-left_x1)/128.;
    T h2 = (right_x2-left_x2)/128.;
    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=left_x1; x1<=right_x1; x1+=h1) {
        //for (T x2=left_x2; x2<=right_x2; x2+=0.03125) {
        for (T x2=left_x2; x2<=right_x2; x2+=h2) {
            T payoff = truncatedOption.payoff(x1,x2);
            T approx = evaluate(basis2d, u, x1, x2);

            maxError = std::max(maxError, fabs(approx-payoff));
            if (fabs(x1)<1. && fabs(x2)<1.) innerMaxError = std::max(innerMaxError, fabs(approx-payoff));
        }
    }

    h1 = 20./128.;
    h2 = 20./128.;    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=-10.; x1<=10.; x1+=h1) {
        for (T x2=-10.; x2<=10.; x2+=h2) {
            T payoff = truncatedOption.payoff(x1,x2);
            T approx = evaluate(basis2d, u, x1, x2);
            plotfile << x1 << " " << x2 << " " << approx << " " << payoff << endl;
        }
        plotfile << endl;
    }
    plotfile.close();
}

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int N,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters)
{
    std::stringstream coeffsfilename;
    coeffsfilename << "cgmye2d_coeffs_endpoint_" << u.size();
    plotScatterCoeff2D(u, basis2d.first, basis2d.second, coeffsfilename.str().c_str());

    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "cgmye2d_basketput_awgm_" << N
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "cgmye2d_sumofputs_awgm_" << N
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);
    std::cerr << "Error is measured over [" << delta*left_x1 << ", " << delta*right_x1 << "], ["
              << delta*left_x2 << ", " << delta*right_x2 << "]" << std::endl;
    T maxerror = 0.;
    T h1 = (delta*right_x1-delta*left_x1)/32.;
    T h2 = (delta*right_x2-delta*left_x2)/32.;
    for (T x1=delta*left_x1; x1<=delta*right_x1; x1+=h1) {
        for (T x2=delta*left_x2; x2<=delta*right_x2; x2+=h2) {
            T S1    = strike*std::exp(x1);
            T S2    = strike*std::exp(x2);
            T exact = option2d.value(processparameters,S1,S2,0);
            //T x1hat = x1+(r-0.5*sigma1*sigma1-0.05698818862)*maturity;
            //T x2hat = x2+(r-0.5*sigma2*sigma2-0.74346053633)*maturity;
            T x1hat = x1+r*maturity;
            T x2hat = x2+r*maturity;
            T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
            T approx =std::exp(-r*maturity)*evaluate(basis2d, u, x1hat, x2hat);
            maxerror = std::max(maxerror, fabs(approx-exact));
        }
    }
    /*
    h1 = 20./128.;
    h2 = 20./128.;
    for (T x1=-10.; x1<=10.; x1+=h1) {
        for (T x2=-10.; x2<=10.; x2+=h2) {
            T S1    = strike*std::exp(x1);
            T S2    = strike*std::exp(x2);
            T exact = option2d.value(processparameters,S1,S2,0);
            //T x1hat = x1+(r-0.5*sigma1*sigma1-0.05698818862)*maturity;
            //T x2hat = x2+(r-0.5*sigma2*sigma2-0.74346053633)*maturity;
            T x1hat = x1+r*maturity;
            T x2hat = x2+r*maturity;
            T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
            T approx =std::exp(-r*maturity)*evaluate(basis2d, u, x1hat, x2hat);
            plotfile << x1 << " " << x2 << " " << exact << " " << approx << " " << payoff << endl;
        }
        plotfile << endl;
    }
    plotfile.close();
    */
    return maxerror;
}

void
getSparseGridIndexSet(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                      Coefficients<Lexicographical,T,Index2D> &u, int j)
{

    int j0_x = basis2d.first.j0;
    int j0_y = basis2d.second.j0;

    std::cerr << "getSparseGridIndexSet called" << std::endl;
    for (int levelsum=0; levelsum<=j; ++levelsum) {
        for (int i1=0; i1<=levelsum; ++i1) {
            for (int i2=0; i2<=levelsum; ++i2) {
                if (i1+i2!=levelsum) continue;

                if ((i1==0) && (i2==0)) {
                    std::cerr << "Scaling function part..." << std::endl;
                    long left_k_x1 = 0, right_k_x1 = 0, left_k_x2 = 0, right_k_x2 = 0;
                    getLeftAndRightTranslationIndices(basis2d.first, left_x1, right_x1, j0_x, XBSpline, left_k_x1, right_k_x1);
                    getLeftAndRightTranslationIndices(basis2d.second,left_x2, right_x2, j0_y, XBSpline, left_k_x2, right_k_x2);
                    std::cerr << "[" << left_x1 << ", " << right_x1 << "], [" << left_x2 << ", " << right_x2 << "]" << std::endl;
                    for (long k1=left_k_x1; k1<=right_k_x1; ++k1) {
                       for (long k2=left_k_x2; k2<=right_k_x2; ++k2) {
                           Index1D row(j0_x,k1,XBSpline);
                           Index1D col(j0_y,k2,XBSpline);
                           u[Index2D(row,col)] = 0.;
                       }
                    }
                }
                else if ((i1!=0) && (i2==0)) {
                    int j1=j0_x+i1-1;
                    long left_k_x1 = 0, right_k_x1 = 0, left_k_x2 = 0, right_k_x2 = 0;
                    getLeftAndRightTranslationIndices(basis2d.first, left_x1, right_x1, j1, XWavelet, left_k_x1, right_k_x1);
                    getLeftAndRightTranslationIndices(basis2d.second,left_x2, right_x2, j0_y, XBSpline, left_k_x2, right_k_x2);
                    for (long k1=left_k_x1; k1<=right_k_x1; ++k1) {
                        for (long k2=left_k_x2; k2<=right_k_x2; ++k2) {
                            Index1D row(j1,k1,XWavelet);
                            Index1D col(j0_y,k2,XBSpline);
                            u[Index2D(row,col)] = 0.;
                        }
                    }
                }
                else if ((i1==0) && (i2!=0)) {
                    int j2=j0_y+i2-1;
                    long left_k_x1 = 0, right_k_x1 = 0, left_k_x2 = 0, right_k_x2 = 0;
                    getLeftAndRightTranslationIndices(basis2d.first, left_x1, right_x1, j0_x, XBSpline, left_k_x1, right_k_x1);
                    getLeftAndRightTranslationIndices(basis2d.second,left_x2, right_x2, j2, XWavelet, left_k_x2, right_k_x2);
                    for (long k1=left_k_x1; k1<=right_k_x1; ++k1) {
                        for (long k2=left_k_x2; k2<=right_k_x2; ++k2) {
                            Index1D row(j0_x,k1,XBSpline);
                            Index1D col(j2,k2,XWavelet);
                            u[Index2D(row,col)] = 0.;
                        }
                    }
                }
                else if ((i1!=0) && (i2!=0)) {
                    int j1=j0_x+i1-1;
                    int j2=j0_y+i2-1;
                    long left_k_x1 = 0, right_k_x1 = 0, left_k_x2 = 0, right_k_x2 = 0;
                    getLeftAndRightTranslationIndices(basis2d.first, left_x1, right_x1, j1, XWavelet, left_k_x1, right_k_x1);
                    getLeftAndRightTranslationIndices(basis2d.second,left_x2, right_x2, j2, XWavelet, left_k_x2, right_k_x2);
                    for (long k1=left_k_x1; k1<=right_k_x1; ++k1) {
                        for (long k2=left_k_x2; k2<=right_k_x2; ++k2) {
                            Index1D row(j1,k1,XWavelet);
                            Index1D col(j2,k2,XWavelet);
                            u[Index2D(row,col)] = 0.;
                        }
                    }
                }
            }
        }
    }

}

void
getLeftAndRightTranslationIndices(const PrimalBasis &basis, T left_x, T right_x,
                                  int j, XType type, long &left_k, long &right_k)
{
    Support<T> supp(left_x,right_x);
    if (type == XBSpline) {
        left_k = 0;
        while (overlap(basis.mra.phi.support(j,left_k),supp)>0) {
            left_k--;
        }
        left_k -= 5;
        right_k = 0;
        while (overlap(basis.mra.phi.support(j,right_k),supp)>0) {
            right_k++;
        }
        right_k += 5;
    }
    else {
        left_k = 0;
        while (overlap(basis.psi.support(j,left_k),supp)>0) {
            left_k--;
        }
        left_k -= 5;
        right_k = 0;
        while (overlap(basis.psi.support(j,right_k),supp)>0) {
            right_k++;
        }
        right_k += 5;
    }
}
