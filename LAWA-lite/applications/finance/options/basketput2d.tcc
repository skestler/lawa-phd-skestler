namespace lawa {

template <typename T>
Option2D<T,BasketPut>::Option2D(void)
    : optionparameters(), N(1000000)

{

}

template <typename T>
Option2D<T,BasketPut>::Option2D(const OptionParameters2D<T,BasketPut> &_optionparameters)
    : optionparameters(_optionparameters), N(1000000)

{
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
}

template <typename T>
T
Option2D<T,BasketPut>::payoff(T S1, T S2) const
{
    return std::max(optionparameters.strike - optionparameters.weight1*S1 - optionparameters.weight2*S2, (T)0.);
}

template <typename T>
T
Option2D<T,BasketPut>::payoff_log(T x1, T x2) const
{
    return optionparameters.strike*std::max(1.-optionparameters.weight1*std::exp(x1)
                                              -optionparameters.weight2*std::exp(x2), (T)0.);
}


template <typename T>
T
Option2D<T,BasketPut>::value(const ProcessParameters2D<T,BlackScholes2D> &processparameters,
                             T S1, T S2, T t)
{
    typedef typename std::map<std::pair<T,T>,T>::const_iterator const_it;
    std::pair<T,T> pairS(S1,S2);

    const_it it = values.find(pairS);

    if (it != values.end()) {
        std::cerr << "Values found." << std::endl;
        return (*it).second;
    }
    else {
        T q11      = processparameters.sigma1*processparameters.sigma1;
        T q22      = processparameters.sigma2*processparameters.sigma2;
        T q12      = processparameters.sigma1*processparameters.sigma2*processparameters.rho;
        T q21      = q12;
        T drift1   = -0.5*q11;
        T drift2   = -0.5*q22;
        T r        = processparameters.r;
        T maturity = optionparameters.maturity;
        ND_Generator2D<T>               nd_generator2d(q11,q12,q21,q22,drift1,drift2);

        T optionPrice_estim = 0.;
        for (int j=1; j<=N; ++j) {
            T x1 = 0., x2 = 0.;
            nd_generator2d(maturity,x1,x2);
            T val = payoff(S1*exp(r*maturity+x1),S2*exp(r*maturity+x2));
            optionPrice_estim += exp(-r*maturity)*val;
        }
        optionPrice_estim /= N;
        values[pairS] = optionPrice_estim;
        std::cout.precision(16);
        std::cout << "Option price: " << optionPrice_estim << std::endl;
        return optionPrice_estim;
    }
}

template <typename T>
T
Option2D<T,BasketPut>::value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters,
                             T S1, T S2, T t)
{
    return 0.;
}

template <typename T>
void
Option2D<T,BasketPut>::setNumberOfMCRuns(int _N) {
    N = _N;
}

}   //namespace lawa
