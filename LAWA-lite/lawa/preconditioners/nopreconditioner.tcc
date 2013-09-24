namespace lawa {

template <typename T, typename Index>
T
NoPreconditioner<T, Index>::operator()(XType /*xtype*/, int /*j*/, int /*k*/) const
{
    return 1.;
}

template <typename T, typename Index>
T
NoPreconditioner<T, Index>::operator()(const Index &/*index*/) const
{
    return 1.;
}

}   // namespace lawa

