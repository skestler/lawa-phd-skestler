#include <fstream>
#include <iostream>
#include <flens/flens.h>
#include <lawa/lawa.h>
#include <applications/finance/fourierpricing/fourierpricing.h>
#include <applications/finance/initialconditions/initialconditions.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>
#include <applications/finance/righthandsides/righthandsides.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

const OptionType1D optiontype = Put;
T strike = 100.;
T maturity = 1.;
T S0 = 100.;
OptionParameters1D<T,Put> optionparameters(strike, maturity, false);

//const ProcessType1D  processtype  = BlackScholes;
//const ProcessType1D  processtype  = CGMY;
const ProcessType1D  processtype  = CGMYe;

/* Reference values for Europ. option from Fang & Oosterlee (2008) and Almendral & Oosterlee (2007)
 * Option parameters strike = 100, maturity = 1, SO = 100
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =1.5   49.790905469 (call)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.8   14.789424 (put)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.1    6.353404 (put)
 */

typedef Basis<T,Primal,Interval,Dijkema>                      Basis1D;
//typedef Basis<T,Orthogonal,Interval,Multi>                      Basis1D;

typedef Integral<Gauss,Basis1D,Basis1D>                         IntegralBasis1DBasis1D;
typedef IntegralF<Gauss,Basis1D>                                IntegralFBasis1D;

typedef IdentityOperator1D<T, Basis1D>                          ScalarProduct1D;

typedef FinanceOperator1D<T, processtype, Basis1D>              FinanceOp;

typedef OptionRHS1D<T, optiontype, processtype, Basis1D>        OptionRhs;

// TimeStepping Methods
typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp, OptionRhs>     ThetaStepScalarProduct1D;
typedef TimeStepping<T, ThetaStepScalarProduct1D>               TimeStepperScalarProduct1D;


template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J, const DenseVectorT &u,
                              const DenseVectorT &u0, T R1, T R2, bool excessToPayoff,
                              T &L2error, T &Linftyerror);

void
getPu0(const Basis1D &basis, DenseVectorT &Pu0,  const Option1D<T,Put> &option, T R1, T R2, int J);

template<typename T, OptionType1D OType, ProcessType1D PType>
void
plotOptionPriceSurface(Option1D<T,OType> &option, const ProcessParameters1D<T,PType> &processparameters, T R1, T R2);

template<typename T, ProcessType1D PType, typename Basis>
void
spyStiffnessMatrix(const Basis &basis, T R1, T R2, const FinanceOp &finance_op, int J, T tol,
                   ProcessParameters1D<T,PType> &processparameters, bool pattern);

int
main(int argc, char *argv[])
{
    cout.precision(12);
    if (argc != 13) {
        cerr << "usage: " << argv[0] << " j0 J R1 R2 excessToPayoff theta timesteps r sigma G M Y" << endl;
        exit(1);
    }
    int d=2, d_=2;
    int j0         = atoi(argv[1]);
    int j_max      = atoi(argv[2]);

    T   R1         = T(atof(argv[3]));
    T   R2         = T(atof(argv[4]));
    int etp        = atoi(argv[5]);

    T theta        = T(atof(argv[6]));
    int timesteps  = atoi(argv[7]);

    T   r          = T(atof(argv[8]));
    T   sigma      = T(atof(argv[9]));
    T   G          = T(atof(argv[10]));
    T   M          = T(atof(argv[11]));
    T   Y          = T(atof(argv[12]));

    bool excessToPayoff   = (etp == 1) ? true : false;

    //ProcessParameters1D<T,BlackScholes>   processparameters(0.04, 0.2);

    //ProcessParameters1D<T,CGMY>           processparameters(r, 1., G, M, Y);
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 1.5 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 0.8 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 0.1 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.04, 1., 2.4, 4.5, 1.8 );

    ProcessParameters1D<T,CGMYe>          processparameters(r, 1., G, M, Y, sigma);
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 2.4, 4.5, 1.8, 0.1 );
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 7.4, 8.5, 1.1, 0.1 );
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 5., 5., 1.5, 0.1 );

    cout << processparameters << endl;

    if (theta < 0.5) {
        cout << "theta should be larger than 0.5!" << endl;
        exit(1);
    }

    T                       timestep = optionparameters.maturity/T(timesteps);
    //Basis1D                 basis(d,d_,j0);
    Basis1D                 basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    int                             order=20;
    Option1D<T,optiontype>  option(optionparameters);

    //plotOptionPriceSurface(option, processparameters, R1, R2);
    //return 0;

    std::stringstream filename;
    if (excessToPayoff) {
        filename << "wavelet_galerkin_option_pricer1d_conv_" << d << "_" << d_
                 << "_R1_" << R1 << "_R2_" << R2 << "_wetp__"
                 << strike << "_" << maturity << "_" << processparameters << ".txt";
    }
    else {
        filename << "wavelet_galerkin_option_pricer1d_conv_" << d << "_" << d_
                 << "_R1_" << R1 << "_R2_" << R2 << "_woetp__"
                 << strike << "_" << maturity << "_" << processparameters << ".txt";
    }
    std::ofstream convfile(filename.str().c_str());

    for (int J=j0; J<=j_max; ++J) {

        Timer time;
        time.start();
        DenseVectorT u(basis.mra.rangeI(J)), u0(basis.mra.rangeI(J));

        if (!excessToPayoff) {
            cout << "Not using excess to payoff!" << endl;
            getPu0(basis, u0, option, R1, R2, J);
        }

        FinanceOp                         finance_op(basis, processparameters, R1, R2, order, J);
        OptionRhs                         rhs(optionparameters, processparameters, basis,
                                              R1, R2, excessToPayoff);
        ThetaStepScalarProduct1D          scheme(theta, basis, finance_op, rhs,
                                                 true, false, 0., 1e-12);

        //spyStiffnessMatrix(basis, R1, R2, finance_op, J, (T)0., processparameters, false);
        //continue;


        TimeStepperScalarProduct1D        timestepmethod(scheme, timestep, timesteps, J);
        u = timestepmethod.solve(u0, false);

        time.stop();
        T comp_time = time.elapsed();

        T L2error = 0.;
        T Linftyerror = 0.;
        ComputeL2ErrorAndPlotSolution(option, processparameters, basis, J, u, u0, R1, R2,
                                      excessToPayoff, L2error, Linftyerror);
        cout      << J << " " << basis.mra.cardI(J) << " " << comp_time << " "
                  << L2error << " " << Linftyerror << endl;
        convfile  << J << " " << basis.mra.cardI(J) << " " << comp_time << " "
                  << L2error << " " << Linftyerror << endl;
    }
    return 0;
}

template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J, const DenseVectorT &u,
                              const DenseVectorT &u0, T R1, T R2, bool excessToPayoff,
                              T &L2error, T &Linftyerror)
{
    std::stringstream filename;
    filename << "tmp2.txt";
    std::ofstream plotFile(filename.str().c_str());


    TruncatedPutOption1D<T,OType> truncatedPutOption;
    truncatedPutOption.setOption(option);
    T h1 = (R1+R2)*pow2i<T>(-7);
    truncatedPutOption.setTruncation(-R1,R2,1,h1);

    T tmp_maturity = option.optionparameters.maturity;
    T tmp_strike   = option.optionparameters.strike;
    T r        = processparameters.r;

    T approxPutS0  = exp(-r*tmp_maturity)*(1./sqrt(R1+R2))*evaluate(basis, J, u, (r*tmp_maturity+R1)/(R1+R2), 0);
    T approxCallS0 = approxPutS0 + S0 - tmp_strike*exp(-processparameters.r*tmp_maturity);
    T exactPutS0   = option.value(processparameters, tmp_strike, 0);
    T exactCallS0  = exactPutS0 + S0 - tmp_strike*exp(-processparameters.r*tmp_maturity);

    std::cout << "Reference value put : " << exactPutS0 << ", approximated value: " << approxPutS0 << std::endl;
    std::cout << "Reference value call: " << exactCallS0 << ", approximated value: " << approxCallS0 << std::endl;

    L2error     = 0.;
    Linftyerror = 0.;
    T h = pow2i<T>(-J-2)*(R1+R2);
    T delta = 0.05;
    //T h = pow2i<T>(-J-2)*(R1+R2);
    //T delta = 1.;
    std::cerr << "Error computation started..." << std::endl;
    for (T x=-delta*R1; x<=delta*R2; x+=h) {
        T approx = (1./sqrt(R1+R2))*evaluate(basis, J, u, (x+R1)/(R1+R2), 0);
        T approx_u0 = (1./sqrt(R1+R2))*evaluate(basis, J, u0, (x+R1)/(R1+R2), 0);
        T exact = 0.;
        T spot = tmp_strike*std::exp(x-r*tmp_maturity);
        exact = std::exp(r*tmp_maturity)*option.value(processparameters, spot, 0);

        if (excessToPayoff) exact -= option.payoff(tmp_strike*exp(x));
        //if (excessToPayoff) approx += option.payoff(tmp_strike*exp(x));

        if ((fabs(x+delta*R1)<1e-12) || (fabs(x-delta*R2) < 1e-12)) {
            L2error += 0.5*h*std::pow(approx-exact,(T)2.);
        }
        else    {
            L2error += h*std::pow(approx-exact,(T)2.);
        }
        Linftyerror = std::max(Linftyerror, fabs(exact-approx));

        plotFile  << x << " " << exact  << " " << approx << " " << option.payoff_log(x) << " "
                  << truncatedPutOption.g_trunc(x) << " " << approx_u0 << endl;

    }
    plotFile.close();
    L2error = sqrt(L2error);
    std::cerr << "...finished." << std::endl;
    //std::cerr << "Please hit enter..." << std::endl;
    //getchar();
}

void
getPu0(const Basis1D &basis, DenseVectorT &Pu0,  const Option1D<T,Put> &option, T R1, T R2, int J)
{
    TruncatedPutOption1D<T,optiontype> truncatedPutOption;
    truncatedPutOption.setOption(option);
    T h = (R1+R2)*pow2i<T>(-7);
    truncatedPutOption.setTruncation(-R1,R2,1,h);

    Function<T> g_fct(truncatedPutOption.g_trunc,truncatedPutOption.singPts);

    InitialCondition1D<Basis1D> initialCondition1D(g_fct,basis,-R1,R2);

    int j0 = basis.j0;
    int N  = basis.mra.cardI(J);

    /*
    DenseVectorT f(basis.mra.rangeI(J));
    int pos = basis.mra.rangeI(J).firstIndex();
    for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        f(pos) = initialCondition1D(XBSpline,j0,k,0);
        ++pos;
    }
    for (int j=j0; j<J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            f(pos) = initialCondition1D(XWavelet,j,k,0);
            ++pos;
        }
    }

    IntegralBasis1DBasis1D integral(basis,basis);
    SparseMatrixT A(N,N);

    int row = 1;
    for (int k_row=basis.mra.rangeI(j0).firstIndex(); k_row<=basis.mra.rangeI(j0).lastIndex(); ++k_row) {
        int col = 1;
        for (int k_col=basis.mra.rangeI(j0).firstIndex(); k_col<=basis.mra.rangeI(j0).lastIndex(); ++k_col) {
            T tmp = integral(j0,k_row,XBSpline,0, j0,k_col,XBSpline,0);
            if (fabs(tmp)>1e-14) A(row,col) = tmp;
            ++col;
        }
        for (int j_col=j0; j_col<J; ++j_col) {
            for (int k_col=basis.rangeJ(j_col).firstIndex(); k_col<=basis.rangeJ(j_col).lastIndex(); ++k_col) {
                T tmp = integral(j0,k_row,XBSpline,0, j_col,k_col,XWavelet,0);
                if (fabs(tmp)>1e-14) A(row,col) = tmp;
                ++col;
            }
        }
        ++row;
    }
    for (int j_row = j0; j_row<J; ++j_row) {
        for (int k_row=basis.rangeJ(j_row).firstIndex(); k_row<=basis.rangeJ(j_row).lastIndex(); ++k_row) {
            int col = 1;
            for (int k_col=basis.mra.rangeI(j0).firstIndex(); k_col<=basis.mra.rangeI(j0).lastIndex(); ++k_col) {
                T tmp = integral(j_row,k_row,XWavelet,0, j0,k_col,XBSpline,0);
                if (fabs(tmp)>1e-14) A(row,col) = tmp;
                ++col;
            }
            for (int j_col=j0; j_col<J; ++j_col) {
                for (int k_col=basis.rangeJ(j_col).firstIndex(); k_col<=basis.rangeJ(j_col).lastIndex(); ++k_col) {
                    T tmp = integral(j_row,k_row,XWavelet,0, j_col,k_col,XWavelet,0);
                    if (fabs(tmp)>1e-14) A(row,col) = tmp;
                    ++col;
                }
            }
            ++row;
        }
    }
    A.finalize();

    int iter = cg(A, Pu0, f, 1e-12,100);
    cerr << "DEBUG: required " << iter << " iterations" << endl;
    */

    int pos = basis.mra.rangeI(J).firstIndex();
    for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        Pu0(pos) = initialCondition1D(XBSpline,j0,k,0);
        ++pos;
    }
    for (int j=j0; j<J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            Pu0(pos) = initialCondition1D(XWavelet,j,k,0);
            ++pos;
        }
    }

}

template<typename T, OptionType1D OType, ProcessType1D PType>
void
plotOptionPriceSurface(Option1D<T,OType> &option, const ProcessParameters1D<T,PType> &processparameters, T R1, T R2)
{
    cout << "plotOptionPriceSurface called." << endl;

    std::stringstream filename;
    filename << "optionpricesurface_" << processparameters << ".txt";
    std::ofstream plotFile(filename.str().c_str());

    Kernel<T,CGMY> kernel(processparameters);

    T tmp_maturity = option.optionparameters.maturity;
    T tmp_strike   = option.optionparameters.strike;
    T r        = processparameters.r;

    T h = 0.2;
    for (T tau=0.1; tau<=tmp_maturity; tau+=0.1) {
        for (T x=-R1; x<=R2; x+=h) {
            plotFile << tau << " " << x << " "
                     << exp(r*tau)*option.value(processparameters,tmp_strike*exp(x-r*tau),tmp_maturity-tau) // kernel.ExpXmOnemX_k_pos+kernel.ExpXmOnemX_k_neg
                     << " " << option.payoff(tmp_strike*exp(x)) << endl;
        }
        plotFile << endl;
    }
}

template<typename T, ProcessType1D PType, typename Basis>
void
spyStiffnessMatrix(const Basis &basis, T R1, T R2, const FinanceOp &finance_op, int J, T tol,
                   ProcessParameters1D<T,PType> &processparameters, bool pattern)
{
    std::stringstream matrixfilename;
    matrixfilename << "A_" << J << "_" << basis.d << "_" << basis.d_ << "_" << tol
                   << "_R1_" << R1 << "_R2_" << R2 << processparameters << ".txt";
    int N =basis.mra.cardI(J);
    SparseMatrixT A(N,N);

    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;

    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
            T val = finance_op(XBSpline, j0, k1, XBSpline, j0, k2);
            if(fabs(val) > tol){
                A(k1-offsetI, k2-offsetI) = pattern ? 1 : val;
            }
            else {
                cout << k1 << ", " << k2 << ": " << val << endl;
            }
        }
    }

    // SF * W
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int j = j0; j <= J-1; ++j){
            for(int k2 = basis.rangeJ(j).firstIndex(); k2 <= basis.rangeJ(j).lastIndex(); ++k2){
                T val = finance_op(XBSpline, j0, k1, XWavelet, j, k2);
                if(fabs(val) > tol){
                    A(k1-offsetI,  basis.mra.cardI(j) + k2 - offsetJ) = pattern ? 1 : val;
                }
            }
        }
    }

      // W * SF
    for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
        for(int j = j0; j <= J-1; ++j){
            for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
                T val = finance_op(XWavelet, j, k1, XBSpline, j0, k2);
                if(fabs(val) > tol){
                    A(basis.mra.cardI(j) + k1 - offsetJ, k2 - offsetI) = pattern ? 1 : val;
                }
            }
        }
    }

        // W * W
    for(int j = j0; j <= J-1; ++j){
        for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
            for(int j_ = j0; j_ <= J-1; ++j_){
                for(int k2 = basis.rangeJ(j_).firstIndex(); k2 <= basis.rangeJ(j_).lastIndex(); ++k2){
                    T val = finance_op(XWavelet, j, k1, XWavelet, j_, k2);
                    if(fabs(val) > tol){
                        A(basis.mra.cardI(j) + k1 - offsetJ, basis.mra.cardI(j_) + k2 - offsetJ) = pattern ? 1 : val;
                    }
                }
            }
        }
    }

    A.finalize();


    spy(A,matrixfilename.str().c_str(),true,(T)0.);
}


