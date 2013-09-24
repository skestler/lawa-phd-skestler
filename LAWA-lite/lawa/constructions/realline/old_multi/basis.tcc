namespace lawa {

template <typename T>
Basis<T,Orthogonal,R,Multi>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _j(j0), psi(*this)
{
}

template <typename T>
int
Basis<T,Orthogonal,R,Multi>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,R,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Orthogonal,R,Multi> &
Basis<T,Orthogonal,R,Multi>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

}   // namespace lawa
