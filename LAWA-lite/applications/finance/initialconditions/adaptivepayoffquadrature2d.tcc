namespace lawa {

template <typename PayoffIntegral>
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::AdaptivePayoffQuadrature2D
(PayoffIntegral &_payoffintegral)
: payoffintegral(_payoffintegral), thresh(1e-14)
{
    rootSolverType = gsl_root_fdfsolver_newton;
    rootSolver = gsl_root_fdfsolver_alloc (rootSolverType);
    FDF.f   = &sing_basketput2d_f;
    FDF.df  = &sing_basketput2d_df;
    FDF.fdf = &sing_basketput2d_fdf;
}

template <typename PayoffIntegral>
typename PayoffIntegral::T
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::integrate(T a1, T b1, T a2, T b2)
{
    T payoff_a1a2 = payoffintegral.payoff(a1,a2);
    T payoff_a1b2 = payoffintegral.payoff(a1,b2);
    T payoff_b1a2 = payoffintegral.payoff(b1,a2);
    T payoff_b1b2 = payoffintegral.payoff(b1,b2);

    if (payoff_a1a2>0. && payoff_a1b2>0. && payoff_b1a2>0. && payoff_b1b2>0.) {
        return this->integrate_smooth(a1, b1, a2, b2);
    }
    if (payoff_a1a2==0 && payoff_a1b2==0 && payoff_b1a2==0 && payoff_b1b2==0) {
        return 0.;
    }

    T y1a2, y1b2, a1y2, b1y2;
    T val = 0.;
    if ( (payoff_a1a2==0 && payoff_a1b2>0) || (payoff_a1a2>0 && payoff_a1b2==0) ) {
        a1y2 = this->find_intersectionpoint_y2_given_y1(a1);
        if ( (payoff_a1a2==0 && payoff_b1a2>0) || (payoff_a1a2>0 && payoff_b1a2==0) )  {
            y1a2 = this->find_intersectionpoint_y1_given_y2(a2);
            val += this->integrate_singular(a1,y1a2,a2,a1y2);
            val += this->integrate_smooth(y1a2,b1,a2,a1y2);
            val += this->integrate_smooth(a1,y1a2,a1y2,b2);
            val += this->integrate_smooth(y1a2,b1,a1y2,b2);
            return val;
        }
        if ( (payoff_a1b2==0 && payoff_b1b2>0) || (payoff_a1b2>0 && payoff_b1b2==0) ) {
            y1b2 = this->find_intersectionpoint_y1_given_y2(b2);
            val += this->integrate_smooth(a1,y1b2,a2,a1y2);
            val += this->integrate_smooth(y1b2,b1,a2,a1y2);
            val += this->integrate_singular(a1,y1b2,a1y2,b2);
            val += this->integrate_smooth(y1b2,b2,a1y2,b2);
            return val;
        }
    }
    if ( (payoff_b1a2==0 && payoff_b1b2>0) || (payoff_b1a2>0 && payoff_b1b2==0) ) {
        b1y2 = this->find_intersectionpoint_y2_given_y1(b1);
        if ( (payoff_a1a2==0 && payoff_b1a2>0) || (payoff_a1a2>0 && payoff_b1a2==0) ) {
            y1a2 = this->find_intersectionpoint_y1_given_y2(a2);
            val += this->integrate_smooth(a1,y1a2,a2,b1y2);
            val += this->integrate_singular(y1a2,b1,a2,b1y2);
            val += this->integrate_smooth(a1,y1a2,b1y2,b2);
            val += this->integrate_smooth(y1a2,b1,b1y2,b2);
            return val;
        }
        if ( (payoff_a1b2==0 && payoff_b1b2>0) && (payoff_a1b2>0 && payoff_b1b2==0) ) {
            y1b2 = this->find_intersectionpoint_y1_given_y2(b2);
            val += this->integrate_smooth(a1,y1b2,a2,b1y2);
            val += this->integrate_smooth(y1b2,b1,a2,b1y2);
            val += this->integrate_smooth(a1,y1b2,b1y2,b2);
            val += this->integrate_singular(y1b2,b2,b1y2,b2);
            return val;
        }
    }



    return 0.;
}

template <typename PayoffIntegral>
typename PayoffIntegral::T
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::integrate_smooth(T a1, T b1, T a2, T b2)
{
    return 0.;
}

template <typename PayoffIntegral>
typename PayoffIntegral::T
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::integrate_singular(T a1, T b1, T a2, T b2)
{
    return 0.;
}

template <typename PayoffIntegral>
bool
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::contains_singularity(T a1, T b1, T a2, T b2)
{
    T weight1 = payoffintegral.option.optionparameters.weight1;
    T weight2 = payoffintegral.option.optionparameters.weight2;
    T u11     = payoffintegral.processparameters.u11;
    T u12     = payoffintegral.processparameters.u12;
    T u21     = payoffintegral.processparameters.u21;
    T u22     = payoffintegral.processparameters.u22;

    T val_a1a2 = weight1*std::exp(u11*a1+u21*a2) + weight2*std::exp(u12*a1+u22*a2) - 1;
    T val_a1b2 = weight1*std::exp(u11*a1+u21*b2) + weight2*std::exp(u12*a1+u22*b2) - 1;
    T val_b1a2 = weight1*std::exp(u11*b1+u21*a2) + weight2*std::exp(u12*b1+u22*a2) - 1;
    T val_b1b2 = weight1*std::exp(u11*b1+u21*b2) + weight2*std::exp(u12*b1+u22*b2) - 1;

    // if the values have different signs, then the singularity enters the integration domain at
    // the corresponding axis
    if (val_a1a2*val_b1a2 < 0) return true;
    if (val_b1a2*val_b1b2 < 0) return true;
    if (val_b1b2*val_a1b2 < 0) return true;
    if (val_a1b2*val_a1a2 < 0) return true;




}

template <typename PayoffIntegral>
typename PayoffIntegral::T
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::find_intersectionpoint_y1_given_y2(T y2)
{
    T weight1 = payoffintegral.option.optionparameters.weight1;
    T weight2 = payoffintegral.option.optionparameters.weight2;
    T u11     = payoffintegral.processparameters.u11;
    T u12     = payoffintegral.processparameters.u12;
    T u21     = payoffintegral.processparameters.u21;
    T u22     = payoffintegral.processparameters.u22;

    double a = weight2*std::exp((1./u22)*y2);
    double b = u21/u22;
    double c = (u21+u22)/u22;
    double d = -weight1;

    std::cerr << "y2 = " << y2 << std::endl;

    struct my_params params = {a, b, c, d};
    std::cerr.precision(16);
    FDF.params = &params;

    int status;
    double x1 = 0.;
    gsl_root_fdfsolver_set (rootSolver, &FDF, x1);

    int iter = 0;
    for (int iter=1; iter<=10; ++iter) {
        status = gsl_root_fdfsolver_iterate (rootSolver);
        x1 = gsl_root_fdfsolver_root (rootSolver);
        double res = sing_basketput2d_f(x1, &params);
        std::cerr <<  iter << " " <<  x1 << " " <<  res << std::endl;
        if (fabs(res)<1e-14) {
            std::cerr << "Newton converged." << std::endl;
            break;
        }
    }
    double x2 = (y2 - u21*x1)/u22;
    double y1 = u11*x1 + u12*x2;

    double res2 = 1. - weight1*std::exp(u11*y1+u21*y2) - weight2*std::exp(u12*y1+u22*y2);
    std::cerr << "(" << y1 << ", " << y2 << "): val = " << res2 << std::endl;
    return y2;
}

template <typename PayoffIntegral>
typename PayoffIntegral::T
AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>::find_intersectionpoint_y2_given_y1(T y1)
{
    T weight1 = payoffintegral.option.optionparameters.weight1;
    T weight2 = payoffintegral.option.optionparameters.weight2;
    T u11     = payoffintegral.processparameters.u11;
    T u12     = payoffintegral.processparameters.u12;
    T u21     = payoffintegral.processparameters.u21;
    T u22     = payoffintegral.processparameters.u22;

    double a = weight2*std::exp((1./u12)*y1);
    double b = u11/u12;
    double c = (u11+u12)/u12;
    double d = -weight1;

    std::cerr << "y1 = " << y1 << std::endl;

    struct my_params params = {a, b, c, d};
    std::cerr.precision(16);
    FDF.params = &params;

    int status;
    double x1 = 0.;
    gsl_root_fdfsolver_set (rootSolver, &FDF, x1);

    int iter = 0;
    for (int iter=1; iter<=10; ++iter) {
        status = gsl_root_fdfsolver_iterate (rootSolver);
        x1 = gsl_root_fdfsolver_root (rootSolver);
        double res = sing_basketput2d_f(x1, &params);
        std::cerr <<  iter << " " <<  x1 << " " <<  res << std::endl;
        if (fabs(res)<1e-14) {
            std::cerr << "Newton converged." << std::endl;
            break;
        }
    }
    double x2 = (y1 - u11*x1)/u12;
    double y2 = u21*x1 + u22*x2;

    double res2 = 1. - weight1*std::exp(u11*y1+u21*y2) - weight2*std::exp(u12*y1+u22*y2);
    std::cerr << "(" << y1 << ", " << y2 << "): val = " << res2 << std::endl;
    return y2;
}

}   // namespace lawa
