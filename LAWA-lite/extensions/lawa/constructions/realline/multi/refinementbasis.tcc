namespace lawa {

template <typename T>
Basis<T,Orthogonal,R,MultiRefinement>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _j(j0), LaplaceOp1D(_d, *this), IdentityOp1D(_d, *this)
{
    if (d > 4) {
        std::cerr << "Basis<T,Orthogonal,R,MultiRefinement> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
    setLevel(_j);
}

template <typename T>
Basis<T,Orthogonal,R,MultiRefinement>::~Basis()
{

}

template <typename T>
int
Basis<T,Orthogonal,R,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,R,MultiRefinement>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Orthogonal,R,MultiRefinement> &
Basis<T,Orthogonal,R,MultiRefinement>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        std::cerr << "BasisFunction<T,Orthogonal,Interval,MultiRefinement> has no wavelet member."
                  << std::endl;
        exit(1);
        return mra.phi;
    }
}



template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Orthogonal,R,MultiRefinement>::
getBSplineNeighborsForBSpline(int j_bspline1, long k_bspline1,
                              const SecondRefinementBasis &secondrefinementbasis,
                              int &j_bspline2, long &k_bspline2_first, long &k_bspline2_last) const
{
    ct_assert(SecondRefinementBasis::Side==Orthogonal and SecondRefinementBasis::Domain==R
              and SecondRefinementBasis::Cons==MultiRefinement);
    //if (flens::IsSame<Basis<T,Orthogonal,Interval,MultiRefinement>, SecondRefinementBasis>::value)
    j_bspline2 = j_bspline1;
    k_bspline2_first = k_bspline1 - d + 1;
    k_bspline2_last  = k_bspline1 + d - 1;

    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,R,MultiRefinement>
::getWaveletNeighborsForBSpline(int j_bspline, long k_bspline, const SecondBasis &secondbasis,
                                int &j_wavelet, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==R
              and SecondBasis::Cons==Multi);

    j_wavelet = j_bspline - secondbasis.psi._addRefinementLevel + 1;
    Support<T> supp = mra.phi.support(j_bspline,k_bspline);
    T a = supp.l1, b = supp.l2;

    long k_tilde = (long)std::floor(pow2i<T>(j_wavelet)*a);
    k_tilde += 1;
    k_wavelet_first  = (k_tilde-2)*secondbasis._numSplines;
    k_wavelet_last   = (k_tilde+2)*secondbasis._numSplines;

    return;
}

template <typename T>
Basis<T,Orthogonal,R,MultiRefinement>::LaplaceOperator1D::
LaplaceOperator1D(int _d, const Basis<T,Orthogonal,R,MultiRefinement> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 1:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D "
            << " will not work d=" << d << std::endl;
            break;
        case 2:
            values1.engine().resize(2,0);
            values1 = 2., -1.;
            break;
        case 3:
            values1.engine().resize(3,0);
            values1 = 1.L, -1.L/3.L, -1.L/6.L;
            break;
        case 4:
            values1.engine().resize(4,0);
            values1 = 6.L/5.L, 0.L, -3.L/8.L, -3.L/40.L;
            values2.engine().resize(4,0);
            values2 = 6.L/5.L, -3.L/8.L, -3.L/8.L, 0.L;
            break;
        default:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
            exit(1);
    }
}

template <typename T>
T
Basis<T,Orthogonal,R,MultiRefinement>::LaplaceOperator1D::
operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2)
{
    assert(j1==j2);
    if (d==2) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  pow2i<T>(2*j1)*values1(k_diff);
    }
    else if (d==3) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        return pow2i<T>(2*j1)*values1(k_diff);
    }
    else if (d==4) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>3) return 0.L;

        int type = k1 > k2 ? refinementbasis.mra.phi._type(k1) : refinementbasis.mra.phi._type(k2);
        if   (type==1)  return pow2i<T>(2*j1)*values1(k_diff);   // Attention: different indexing when compared to interval
        else            return pow2i<T>(2*j1)*values2(k_diff);
    }
    else {
        std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

template <typename T>
Basis<T,Orthogonal,R,MultiRefinement>::IdentityOperator1D::
IdentityOperator1D(int _d, const Basis<T,Orthogonal,R,MultiRefinement> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 1:
            values1.engine().resize(2,0);
            values1 = 0.L, 0.L;
            break;
        case 2:
            values1.engine().resize(2,0);
            values1 = 2.L/3.L, 1.L/6.L;
            break;
        case 3:
            values1.engine().resize(3,0);
            values1 = 11.L/20.L, 13.L/60.L, 1.L/120.L;
            break;
        case 4:
            values1.engine().resize(4,0);
            values1 = 8.L/35.L, 1.L/7.L, 9.L/560.L, 1.L/560.L;
            values2.engine().resize(4,0);
            values2 = 8.L/35.L, 53.L/560.L, 9.L/560.L, 0.L;
            break;
        default:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
            exit(1);
    }
}

template <typename T>
T
Basis<T,Orthogonal,R,MultiRefinement>::IdentityOperator1D::
operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2)
{
    assert(j1==j2);
    if (d==2) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  values1(k_diff);
    }
    else if (d==3) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        return values1(k_diff);
    }
    else if (d==4) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>3) return 0.L;

        int type = k1 > k2 ? refinementbasis.mra.phi._type(k1) : refinementbasis.mra.phi._type(k2);
        if (type==1)  return values1(k_diff);  // Attention: different indexing when compared to interval
        else          return values2(k_diff);
    }
    else {
        std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

} // namespace lawa

