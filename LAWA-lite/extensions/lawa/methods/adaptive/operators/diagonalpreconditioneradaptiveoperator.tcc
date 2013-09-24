namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator>
DiagonalPreconditionerAdaptiveOperator<T,Index,AdaptiveOperator>::
DiagonalPreconditionerAdaptiveOperator(AdaptiveOperator &_A)
    : A(_A)
{

}

template <typename T, typename Index, typename AdaptiveOperator>
T
DiagonalPreconditionerAdaptiveOperator<T,Index,AdaptiveOperator>::operator()(const Index &index)
{
    return A.prec(index);
}

}   // namespace lawa
