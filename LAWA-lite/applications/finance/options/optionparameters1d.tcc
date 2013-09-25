namespace lawa {

template <typename T>
OptionParameters1D<T,Put>::OptionParameters1D(void)
: strike((T)0.), maturity((T)0.), earlyExercise(false)
{

}

template <typename T>
OptionParameters1D<T,Put>::OptionParameters1D(T _strike, T _maturity, bool _earlyExercise)
: strike(_strike), maturity(_maturity), earlyExercise(_earlyExercise)
{

}

}   // namespace lawa
