namespace lawa {

template <typename T, typename Index, typename Basis>
NoCompression<T,Index,Basis>::NoCompression(const Basis &_basis)
    : basis(_basis)
{
}

template <typename T, typename Index, typename Basis>
void
NoCompression<T,Index,Basis>::setParameters(const IndexSet<Index> &/*LambdaRow*/)
{

}

template <typename T, typename Index, typename Basis>
IndexSet<Index>
NoCompression<T,Index,Basis>::SparsityPattern(const Index &/*lambda_col*/, const IndexSet<Index> &LambdaRow)
{
    return LambdaRow;
}

}    //namespace lawa

