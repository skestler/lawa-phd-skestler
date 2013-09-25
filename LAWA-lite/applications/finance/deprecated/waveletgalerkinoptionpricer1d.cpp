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

const ProcessType1D  processtype  = BlackScholes;
ProcessParameters1D<T,BlackScholes>   processparameters(0.04, 0.2);

/* Reference values for Europ. option from Fang & Oosterlee (2008) and Almendral & Oosterlee (2007)
 * Option parameters strike = 100, maturity = 1, SO = 100
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =1.5   49.790905469 (call)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.8   14.789424 (put)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.1    6.353404 (put)
 */

//const ProcessType1D  processtype  = CGMY;
//ProcessParameters1D<T,CGMY>             processparameters(0.1, 1., 5., 5., 1.5 );
//ProcessParameters1D<T,CGMY>             processparameters(0.1, 1., 5., 5., 0.8 );
//ProcessParameters1D<T,CGMY>             processparameters(0.1, 1., 5., 5., 0.1 );
//ProcessParameters1D<T,CGMY>             processparameters(0.04, 1., 2.4, 4.5, 1.8 );

//const ProcessType1D  processtype  = CGMYe;
//ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 2.4, 4.5, 1.8, 0.1 );
//ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 5., 5., 1.5, 0.1 );


typedef Basis<T,Primal,Interval,Dijkema>        Basis1D;

typedef IdentityOperator1D<T, Basis1D>                          ScalarProduct1D;
typedef WeightedIdentityOperator1D<T, Basis1D>                  WeightedScalarProduct1D;
typedef DiagonalMatrixPreconditioner1D<T, Basis1D,
                                      WeightedScalarProduct1D>  WeightedL2Preconditioner1D;

typedef FinanceOperator1D<T, processtype, Basis1D>              FinanceOp;

typedef OptionRHS1D<T, optiontype, processtype, Basis1D>        OptionRhs;

typedef PayoffInitialCondition1D<optiontype, Basis1D>           PayoffInitCond;

// TimeStepping Methods
typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp, OptionRhs>     ThetaStepScalarProduct1D;
typedef TimeStepping<T, ThetaStepScalarProduct1D>               TimeStepperScalarProduct1D;

typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp,
               HomogeneousRHS<T>, WeightedScalarProduct1D >     ThetaStepWeightedScalarProduct1DHomRHS;
typedef TimeStepping<T,
               ThetaStepWeightedScalarProduct1DHomRHS>          TimeStepperWeightedScalarProduct1DHomRHS;

const T                         eta=2.;
ExponentialWeightFunction1D<T>  exponentialweight;

template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J,
                              const DenseVectorT &u0, const DenseVectorT &u,
                              T &L2error, T &Linftyerror, T R1, T R2, int use_excess_to_payoff);

int
main(int argc, char *argv[])
{
    cout.precision(12);
    if (argc != 6) {
        cerr << "usage: " << argv[0] << " j0 J theta timesteps use_excess_to_payoff" << endl;
        exit(1);
    }
    int d=2, d_=2;
    int j0                   = atoi(argv[1]);
    int j_max                = atoi(argv[2]);
    T theta                  = T(atof(argv[3]));
    int timesteps            = atoi(argv[4]);
    int use_excess_to_payoff = atoi(argv[5]);

    if (theta < 0.5) {
        cout << "theta should be larger than 0.5!" << endl;
        exit(1);
    }

    exponentialweight.setEta(eta);
    T                       timestep = optionparameters.maturity/T(timesteps);
    T                       R1=6.;
    T                       R2=6.;
    Basis1D                 basis(d,d_,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    int                             order=20;
    Option1D<T,optiontype>  option(optionparameters);


    cout << "Reference option value put : " << option.value(processparameters, S0, 0) << endl;
    cout << "Reference option value call: " << option.value(processparameters, S0, 0) + S0 - strike*exp(-processparameters.r*maturity) << endl;

    std::stringstream filename;
    filename << "wavelet_galerkin_option_pricer1d_conv_" << R1 << "_" << R2 << "_" << strike << "_"
                 << maturity << "_" << processparameters.sigma << ".txt";
    //filename << "wavelet_galerkin_option_pricer1d_conv_" << R1 << "_" << R2 << "_" << strike << "_"
    //         << maturity << "_" << processparameters.k_C << "_" << processparameters.k_G << "_"
    //         << processparameters.k_M << "_" << processparameters.k_Y << ".txt";
    std::ofstream convfile(filename.str().c_str());

    for (int J=j0; J<=j_max; ++J) {
        Timer time;
        time.start();
        DenseVectorT u(basis.mra.rangeI(J)), u0(basis.mra.rangeI(J));
        if (use_excess_to_payoff==1) {
            cout << "Initializing operator..." << endl;
            FinanceOp                         finance_op(basis, processparameters, 0., R1, R2, order, J);
            cout << "... finished." << endl;
            cout << "Initializing rhs..." << endl;
            OptionRhs                         rhs(optionparameters, processparameters, basis,
                                                  R1, R2);
            cout << "... finished." << endl;
            //ScalarProduct1D                   scalarproduct(basis);
            cout << "Initializing time stepper..." << endl;
            ThetaStepScalarProduct1D          scheme(theta, basis, finance_op, rhs,
                                                     true, false, 0., 1e-12);
            cout << "... finished." << endl;
            cout << "Initializing timestepping scheme..." << endl;
            TimeStepperScalarProduct1D        timestepmethod(scheme, timestep, timesteps, J);
            cout << "... finished." << endl;
            cout << "Time step solver started..." << endl;
            u = timestepmethod.solve(u0, false);
            cout << "... finished." << endl;
        }
        else {
            Function<T>                    weightFct(exponentialweight.weight,
                                                     exponentialweight.singularPoints);
            WeightedScalarProduct1D        wL2scalarproduct(basis, weightFct, order, -R1, R2);
            WeightedL2Preconditioner1D     wpreconditioner(wL2scalarproduct);


            PayoffInitCond                 payoff_initcond(option,basis,weightFct,eta,-R1,R2);
            Assembler1D<T, Basis1D>        assembler(basis);
            SparseMatrixT   M      =       assembler.assembleStiffnessMatrix(wL2scalarproduct, J);
            DiagonalMatrixT P      =       assembler.assemblePreconditioner(wpreconditioner, J);
            DenseVectorT    rhs_u0 =       assembler.assembleRHS(payoff_initcond, J);

            cout << pcg(P, M, u0, rhs_u0) << " pcg iterations required." << endl;

            FinanceOp                 finance_op(basis, processparameters, eta, R1, R2, order, J);
            HomogeneousRHS<T>         homogeneousrhs;
            ThetaStepWeightedScalarProduct1DHomRHS scheme(theta, basis, finance_op, homogeneousrhs,
                                                          wL2scalarproduct, true, false, 0., 1e-12);

            TimeStepperWeightedScalarProduct1DHomRHS timestepmethod(scheme, timestep, timesteps, J);
            u = timestepmethod.solve(u0, false);

        }
        time.stop();
        T comp_time = time.elapsed();

        T L2error = 0.;
        T Linftyerror = 0.;
        ComputeL2ErrorAndPlotSolution(option, processparameters, basis, J, u0, u,
                                      L2error, Linftyerror, R1, R2, use_excess_to_payoff);
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
                              const Basis &basis, int J,
                              const DenseVectorT &u0, const DenseVectorT &u,
                              T &L2error, T &Linftyerror, T R1, T R2, int use_excess_to_payoff)
{
    std::stringstream filename;
    filename << "tmp.txt";
    std::ofstream plotFile(filename.str().c_str());
    std::stringstream filename2;
    filename2 << "tmp2.txt";
    std::ofstream plotFile2(filename2.str().c_str());
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
    T delta = 0.05;//pow2i<T>(0);
    std::cerr << "Error computation started..." << std::endl;
    for (T x=-delta*R1; x<=delta*R2; x+=h) {
        T P_u0   = (1./sqrt(R1+R2))*evaluate(basis, J, u0, (x+R1)/(R1+R2), 0);
        T approx = (1./sqrt(R1+R2))*evaluate(basis, J, u, (x+R1)/(R1+R2), 0);
        T exact = 0.;
        T spot = tmp_strike*std::exp(x-r*tmp_maturity);
        exact = std::exp(r*tmp_maturity)*option.value(processparameters, spot, 0);

        if (use_excess_to_payoff==1) exact -= option.payoff(tmp_strike*exp(x));

        if ((fabs(x+delta*R1)<1e-12) || (fabs(x-delta*R2) < 1e-12)) {
            L2error += 0.5*h*std::pow(approx-exact,2.);
        }
        else    {
            L2error += h*std::pow(approx-exact,2.);
        }
        Linftyerror = std::max(Linftyerror, fabs(exact-approx));

        plotFile  << x << " " << approx << " " << exact << " "
                  << P_u0 << " " << option.payoff_log(x) << endl;

        T weight = exponentialweight.weight(x);
        plotFile2 << x << " " << weight*approx << " " << weight*exact << " "
                  << weight*P_u0 << " " << weight*option.payoff_log(x) << endl;
    }
    plotFile.close();
    plotFile2.close();
    L2error = sqrt(L2error);
    std::cerr << "...finished." << std::endl;
}



