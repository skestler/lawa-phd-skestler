#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC 1

#include <cassert>
#include <lawa/math/linspace.h>
#include <lawa/math/pow2.h>
#include <lawa/constructions/interval/multi/_constant_evaluator.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _j(j0), LaplaceOp1D(_d, *this), IdentityOp1D(_d, *this)
{
    assert(d>=1);
    setLevel(_j);
    if (d>1) {
        this->enforceBoundaryCondition<DirichletBC>();
    }
}

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::~Basis()
{

}

template <typename T>
int
Basis<T,Orthogonal,Interval,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,MultiRefinement>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Orthogonal,Interval,MultiRefinement>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Orthogonal,Interval,MultiRefinement> &
Basis<T,Orthogonal,Interval,MultiRefinement>::generator(XType xtype) const
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
Basis<T,Orthogonal,Interval,MultiRefinement>::
getBSplineNeighborsForBSpline(int j_bspline1, long k_bspline1,
                              const SecondRefinementBasis &secondrefinementbasis,
                              int &j_bspline2, long &k_bspline2_first, long &k_bspline2_last) const
{
    ct_assert(SecondRefinementBasis::Side==Orthogonal and SecondRefinementBasis::Domain==Interval
              and SecondRefinementBasis::Cons==MultiRefinement);
    //if (flens::IsSame<Basis<T,Orthogonal,Interval,MultiRefinement>, SecondRefinementBasis>::value)
    j_bspline2 = j_bspline1;
    Support<T> supp = mra.phi.support(j_bspline1,k_bspline1);
    T a = supp.l1, b = supp.l2;
    if (a==0.L) {
        k_bspline2_first = 0;       // In this case, range always starts with 0
        k_bspline2_last  = std::min(k_bspline1 + d, (long)mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_bspline2_first = std::max(k_bspline1 - d + 1, 0L);
        k_bspline2_last  = std::min(k_bspline1 + d - 1, (long)mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    k_bspline2_first = std::max(0L,k_bspline1 - d);
    k_bspline2_last  = (long)mra.rangeI(j_bspline2).lastIndex();

    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,MultiRefinement>
::getWaveletNeighborsForBSpline(int j_bspline, long k_bspline, const SecondBasis &secondbasis,
                                int &j_wavelet, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);

    j_wavelet = j_bspline - secondbasis._addRefinementLevel + 1;
    Support<T> supp = mra.phi.support(j_bspline,k_bspline);
    T a = supp.l1, b = supp.l2;

    if (a==0.L) {
        k_wavelet_first = 1;
        k_wavelet_last =  secondbasis._numLeftParts + secondbasis._numInnerParts;
        k_wavelet_last  = std::min(k_wavelet_last, (long)secondbasis.rangeJR(j_wavelet).lastIndex());
        return;
    }
    if (0<a && b<1.L) {
        long k_tilde = (long)std::floor(pow2i<T>(j_wavelet)*a);
        k_tilde += 1;
        T tmp = k_tilde*pow2i<T>(-j_wavelet);
        if (a<tmp && tmp<b) {
            k_wavelet_first  = (k_tilde-2)*secondbasis._numInnerParts;
            k_wavelet_last   = (k_tilde+2)*secondbasis._numInnerParts;
        }
        else {
            k_wavelet_first  = (k_tilde-2)*secondbasis._numInnerParts;
            k_wavelet_last   = (k_tilde+1)*secondbasis._numInnerParts;
        }
        k_wavelet_first = std::max((long)secondbasis.rangeJL(j_wavelet).firstIndex(),k_wavelet_first);
        k_wavelet_last  = std::min((long)secondbasis.rangeJR(j_wavelet).lastIndex(), k_wavelet_last);
        return;
    }
    k_wavelet_last   = secondbasis.rangeJ(j_wavelet).lastIndex();
    k_wavelet_first  = k_wavelet_last - (secondbasis._numRightParts + secondbasis._numInnerParts) + 1;
    k_wavelet_first  = std::max(1L, k_wavelet_first);
    return;
}

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D::
LaplaceOperator1D(int _d, const Basis<T,Orthogonal,Interval,MultiRefinement> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 1:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D "
            << " will not work d=" << d << std::endl;
            break;
        case 2:
            inner_values1.engine().resize(2,0);
            inner_values1 = 2., -1.;
            break;
        case 3:
            outer_values.engine().resize(3,0);
            outer_values = 4.L/3.L, -1.L/6.L, -1.L/6.L;
            inner_values1.engine().resize(3,0);
            inner_values1 = 1.L, -1.L/3.L, -1.L/6.L;
            break;
        case 4:
            outer_values.engine().resize(4,0);
            outer_values = 6.L/5.L, 0.L, -3.L/10.L, 0.L;
            inner_values1.engine().resize(4,0);
            inner_values1 = 6.L/5.L, 0.L, -3.L/8.L, -3.L/40.L;
            inner_values2.engine().resize(4,0);
            inner_values2 = 6.L/5.L, -3.L/8.L, -3.L/8.L, 0.L;
            break;
        default:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
            exit(1);
    }
}

template <typename T>
T
Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D::
operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2)
{
    assert(j1==j2);
    if (d==2) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  pow2i<T>(2*j1)*inner_values1(k_diff);
    }
    else if (d==3) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        int firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        int lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return pow2i<T>(2*j1)*outer_values(k_diff);
        }
        else {
            return pow2i<T>(2*j1)*inner_values1(k_diff);
        }
    }
    else if (d==4) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>3) return 0.L;

        int firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        int lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return pow2i<T>(2*j1)*outer_values(k_diff);
        }
        else {
            if (std::max(k1,k2)%2==0)  return pow2i<T>(2*j1)*inner_values1(k_diff);
            else                       return pow2i<T>(2*j1)*inner_values2(k_diff);
        }
    }
    else {
        std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::LaplaceOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D::
IdentityOperator1D(int _d, const Basis<T,Orthogonal,Interval,MultiRefinement> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 1:
            inner_values1.engine().resize(2,0);
            inner_values1 = 0.L, 0.L;
            break;
        case 2:
            inner_values1.engine().resize(2,0);
            inner_values1 = 2.L/3.L, 1.L/6.L;
            break;
        case 3:
            outer_values.engine().resize(3,0);
            outer_values = 1.L/3.L, 5.L/24.L, 1.L/120.L;
            inner_values1.engine().resize(3,0);
            inner_values1 = 11.L/20.L, 13.L/60.L, 1.L/120.L;
            break;
        case 4:
            outer_values.engine().resize(4,0);
            outer_values = 3.L/35.L, 11.L/140.L, 1.L/70.L, 0.L;
            inner_values1.engine().resize(4,0);
            inner_values1 = 8.L/35.L, 1.L/7.L, 9.L/560.L, 1.L/560.L;
            inner_values2.engine().resize(4,0);
            inner_values2 = 8.L/35.L, 53.L/560.L, 9.L/560.L, 0.L;
            break;
        default:
            std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
            exit(1);
    }
}

template <typename T>
T
Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D::
operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2)
{
    assert(j1==j2);
    if (d==2) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  inner_values1(k_diff);
    }
    else if (d==3) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        int firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        int lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return outer_values(k_diff);
        }
        else {
            return inner_values1(k_diff);
        }
    }
    else if (d==4) {
        long k_diff = std::abs(k1 - k2);
        if (k_diff>3) return 0.L;

        int firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        int lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return outer_values(k_diff);
        }
        else {
            if (std::max(k1,k2)%2==0)  return inner_values1(k_diff);
            else                       return inner_values2(k_diff);
        }
    }
    else {
        std::cerr << "Basis<T,Orthogonal,Interval,MultiRefinement>::IdentityOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC
