namespace lawa {

template <typename T>
MRA<T,Primal,R,SparseMulti>::MRA(int _d, int j)
    : d(_d), j0(j), _j(j0), phi(d)
{

}

template <typename T>
int
MRA<T,Primal,R,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,R,SparseMulti>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

}   // namespace lawa
