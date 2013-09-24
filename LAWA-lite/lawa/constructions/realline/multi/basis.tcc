namespace lawa {

template <typename T>
Basis<T,Orthogonal,R,Multi>::Basis(int _d, int j)
    : mra(_d, j), d(_d), d_(_d), j0(mra.j0), _j(j0), psi(*this), refinementbasis(d,j),
     _numSplines(psi._numSplines)//, _addRefinementLevel(psi._addRefinementLevel),
     //_shiftFactor(psi._shiftFactor)
{
    if (d > 4) {
        std::cerr << "Basis<T,Orthogonal,R,Multi> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
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

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getScalingNeighborsForScaling(int j_scaling1, long k_scaling1,
                                const SecondBasis &secondbasis,
                                int &j_scaling2, long &k_scaling_first, long &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_scaling2 = j_scaling1;
    k_scaling_first = k_scaling1 - 2* mra.phi._numSplines;
    k_scaling_last  = k_scaling1 + 2* mra.phi._numSplines;
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getWaveletNeighborsForScaling(int j_scaling, long k_scaling,
                                const SecondBasis &secondbasis,
                                int &j_wavelet, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_wavelet = j_scaling;
    k_wavelet_first = k_scaling - 2*psi._numSplines;
    k_wavelet_last  = k_scaling + 2*psi._numSplines;
    return;
}

template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Orthogonal,R,Multi>
::getBSplineNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                const SecondRefinementBasis &secondrefinementbasis,
                                int &j_bspline, long &k_bspline_first, long &k_bspline_last) const
{
    ct_assert(SecondRefinementBasis::Side==Orthogonal and SecondRefinementBasis::Domain==R
              and SecondRefinementBasis::Cons==MultiRefinement);
    j_bspline = j_wavelet + psi._addRefinementLevel - 1;
    Support<T> supp = psi.support(j_wavelet,k_wavelet);
    Support<T> refinementspline_refsupp = refinementbasis.mra.phi.support(0,0);
    int factor = d == 4 ? 2 : 1;
    k_bspline_first = factor*(pow2i<T>(j_bspline)*supp.l1 - refinementspline_refsupp.l2);
    k_bspline_last  = factor*(pow2i<T>(j_bspline)*supp.l2 - refinementspline_refsupp.l1);
}


template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getScalingNeighborsForWavelet(int j_wavelet, long k_wavelet,
                                const SecondBasis &secondbasis,
                                int &j_scaling, long &k_scaling_first, long &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_scaling = j_wavelet;
    long k_tilde = 0;
    if (k_wavelet >=0)  k_tilde = (k_wavelet / psi._numSplines) * mra.phi._numSplines;
    else                k_tilde = (k_wavelet / psi._numSplines) * mra.phi._numSplines;
    k_scaling_first = k_tilde - 2* mra.phi._numSplines - 2;
    k_scaling_last  = k_tilde + 2* mra.phi._numSplines + 2;
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1, const SecondBasis &secondbasis,
                                int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1;
    long k_tilde = k_wavelet1;
    k_wavelet_first = k_tilde - 2*psi._numSplines;
    k_wavelet_last  = k_tilde + 2*psi._numSplines;
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getLowerWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                     const SecondBasis &secondbasis,
                                     int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1-1;
    long k_tilde = k_wavelet1 / 2;
    k_wavelet_first = k_tilde - 2*psi._numSplines;
    k_wavelet_last  = k_tilde + 2*psi._numSplines;
    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,Multi>
::getHigherWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                     const SecondBasis &secondbasis,
                                     int &j_wavelet2, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);
    j_wavelet2 = j_wavelet1+1;
    long k_tilde = 2*k_wavelet1;
    k_wavelet_first = k_tilde - 4*psi._numSplines;
    k_wavelet_last  = k_tilde + 3*psi._numSplines;
    return;
}

}   // namespace lawa
