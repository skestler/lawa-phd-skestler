namespace lawa {

template <typename T>
Option1D<T,Put>::Option1D()
    : optionparameters()

{

}

template <typename T>
Option1D<T,Put>::Option1D(const OptionParameters1D<T,Put> &_optionparameters)
    : optionparameters(_optionparameters)

{
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
    singularPoints.engine().resize(1);
    singularPoints(1) = 0.;
}

template <typename T>
T
Option1D<T,Put>::payoff(T S) const
{
    return std::max( optionparameters.strike - S, (T)0.);
}

template <typename T>
T
Option1D<T,Put>::payoff_log(T x) const
{
    return optionparameters.strike*std::max(1.-std::exp(x), (T)0.);
}


template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,BlackScholes> &processparameters, T S, T t) const
{
    assert(optionparameters.earlyExercise==false);
    static boost::math::normal norm;
    T r        = processparameters.r;
    T sigma    = processparameters.sigma;
    T strike   = optionparameters.strike;
    T maturity = optionparameters.maturity;
    T d_1 = ( log(S/strike) + ( r + sigma*sigma/2.0 ) * (maturity-t) )
                 / ( sigma * sqrt(maturity-t) );
    T d_2 = d_1 - sigma * sqrt(maturity-t);
    return strike*exp(-r*(maturity-t))*cdf(norm, -d_2) - S * cdf(norm, -d_1);
}

template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,CGMY> &processparameters, T S, T t)
{
    typedef typename std::map<T,T>::const_iterator const_it;

    assert(optionparameters.earlyExercise==false);
    assert(processparameters.k_Y!=0);
    assert(processparameters.k_Y!=1);
    CharacteristicFunction1D<T,CGMY> charfunc(processparameters);

    FourierPricer1D<T, CGMY> fp(charfunc, S, optionparameters.maturity,
                                std::max(optionparameters.strike-10.,(T)0.),
                                optionparameters.strike+10.);

    if (t==0) {
        const_it it_option = values.find(S);
        if (it_option!=values.end()) return (*it_option).second;
        else {
            fp.solve(10000,20);
            T val = fp(optionparameters.strike);
            values[S] = val;
            return val;
        }
    }
    else {
        //fp.solve(10000,17);
        fp.solve(10000,20);
        //fp.solve(10000,23);
    }
    return fp(optionparameters.strike);
}

template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,CGMYe> &processparameters, T S, T t)
{
    typedef typename std::map<T,T>::const_iterator const_it;

    assert(optionparameters.earlyExercise==false);
    assert(processparameters.k_Y!=0);
    assert(processparameters.k_Y!=1);
    CharacteristicFunction1D<T,CGMYe> charfunc(processparameters);

    FourierPricer1D<T, CGMYe> fp(charfunc, S, optionparameters.maturity-t,
                                 std::max(optionparameters.strike-10.,(T)0.1),
                                 optionparameters.strike+10.);

    if (t==0) {
        const_it it_option = values.find(S);
        if (it_option!=values.end()) return (*it_option).second;
        else {
            std::cerr << "   -> Option1D<T,Put>: No high precision is used for computation of reference"
                      << " values!" << std::endl;
            fp.solve(10000,16);
            //fp.solve(10000,20);
            T val = fp(optionparameters.strike);
            values[S] = val;
            return val;
        }
    }
    else {
        //fp.solve(40000,17);
        fp.solve(10000,20);
        //fp.solve(10000,20);
        //fp.solve(10000,16);
    }
    return fp(optionparameters.strike);
}

}   // namespace lawa
