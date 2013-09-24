namespace lawa {

template <typename T>
MRA<T,Orthogonal,R,MultiRefinement>::MRA(int _d, int j)
    : d(_d), j0(j), _j(j0), phi(d)
{
    if (d > 4) {
        std::cerr << "MRA<T,Orthogonal,R,MultiRefinement> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
}

template <typename T>
int
MRA<T,Orthogonal,R,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,R,MultiRefinement>::setLevel(int j) const
{
    assert(d==1 || j>=j0);
    _j = j;
}

} // namespace lawa

