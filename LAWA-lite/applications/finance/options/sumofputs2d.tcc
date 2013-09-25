namespace lawa {

template <typename T>
Option2D<T,SumOfPuts>::Option2D(void)
    : optionparameters(), optionparameters1(), optionparameters2(),
      option1(), option2()

{

}

template <typename T>
Option2D<T,SumOfPuts>::Option2D(const OptionParameters2D<T,SumOfPuts> &_optionparameters)
    : optionparameters(_optionparameters),
      optionparameters1(_optionparameters.strike1, _optionparameters.maturity, _optionparameters.earlyExercise),
      optionparameters2(_optionparameters.strike2, _optionparameters.maturity, _optionparameters.earlyExercise),
      option1(optionparameters1), option2(optionparameters2)

{
    assert(optionparameters.strike1>0.);
    assert(optionparameters.strike2>0.);
    assert(optionparameters.maturity>=0.);
}

template <typename T>
T
Option2D<T,SumOfPuts>::payoff(T S1, T S2) const
{
    return    optionparameters.weight1*std::max(optionparameters.strike1 - S1,(T)0.)
            + optionparameters.weight2*std::max(optionparameters.strike2 - S2,(T)0.);
}

template <typename T>
T
Option2D<T,SumOfPuts>::payoff_log(T x1, T x2) const
{
    return   optionparameters.weight1 * optionparameters.strike1*std::max(1.-std::exp(x1),(T)0.)
           + optionparameters.weight2 * optionparameters.strike2*std::max(1.-std::exp(x2),(T)0.);
}

template <typename T>
T
Option2D<T,SumOfPuts>::value(const ProcessParameters2D<T,BlackScholes2D> &processparameters,
                             T S1, T S2, T t) const
{
    ProcessParameters1D<T,BlackScholes> processparameters1(processparameters.r, processparameters.sigma1);
    ProcessParameters1D<T,BlackScholes> processparameters2(processparameters.r, processparameters.sigma2);

    return    optionparameters.weight1*option1.value(processparameters1,S1,t)
            + optionparameters.weight2*option2.value(processparameters2,S2,t);
}

template <typename T>
T
Option2D<T,SumOfPuts>::value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters,
                             T S1, T S2, T t)
{
    return    optionparameters.weight1*option1.value(processparameters.proc_param1,S1,t)
            + optionparameters.weight2*option2.value(processparameters.proc_param2,S2,t);
}

}   //namespace lawa
