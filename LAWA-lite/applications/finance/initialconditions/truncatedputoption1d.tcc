namespace lawa {

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::left;

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::right;

template <typename T, OptionType1D OType>
int
TruncatedPutOption1D<T, OType>::type;

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::truncWidth;

template <typename T, OptionType1D OType>
flens::DenseVector<Array<T> >
TruncatedPutOption1D<T, OType>::singPts;

template <typename T, OptionType1D OType>
Option1D<T,OType>
TruncatedPutOption1D<T, OType>::option;

template <typename T, OptionType1D OType>
void
TruncatedPutOption1D<T, OType>::setOption(const Option1D<T,Put> &_option)
{
    option.optionparameters = _option.optionparameters;
    option.singularPoints   = _option.singularPoints;
}

template <typename T, OptionType1D OType>
void
TruncatedPutOption1D<T, OType>::setTruncation(T _left, T _right, int _type, T _truncWidth)
{
    left       = _left;
    right      = _right;
    type       = _type;
    truncWidth = _truncWidth;
    int p = option.singularPoints.engine().length();
    singPts.engine().resize(p+1);

    singPts(1) = left+truncWidth;
    for (int i=1; i<=p; ++i) {
        singPts(1+i) = option.singularPoints(i);
    }
}

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::g_trunc(T x)
{
    if (x>=left+truncWidth) {
        return option.payoff_log(x);
    }
    else  {
        T b = left+truncWidth;
        T a = left;
        T val = option.payoff_log(b);
        return val*(x-a)/(b-a);
    }
}

}   // namespace lawa
