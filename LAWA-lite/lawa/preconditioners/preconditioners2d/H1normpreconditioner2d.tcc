#include <cmath>

namespace lawa {

template <typename T, typename Basis2D>
H1NormPreconditioner2D<T,Basis2D>::H1NormPreconditioner2D(const Basis2D &basis)
    : _integral_x(basis.first, basis.first), _integral_y(basis.second, basis.second),
      refval_id_bspline(0), refval_dd_bspline(0), refval_id_wavelet(0), refval_dd_wavelet(0)
{
    if (    (flens::IsSame<Basis_x, Basis<T,Primal,R,CDF> >::value)
         && (flens::IsSame<Basis_y, Basis<T,Primal,R,CDF> >::value) ){
        refval_id_bspline = _integral_x(0,0,XBSpline,0,0,0,XBSpline,0);
        refval_dd_bspline = _integral_x(0,0,XBSpline,1,0,0,XBSpline,1);
        refval_id_wavelet = _integral_x(0,0,XWavelet,0,0,0,XWavelet,0);
        refval_dd_wavelet = _integral_x(0,0,XWavelet,1,0,0,XWavelet,1);
    }
}

template <typename T, typename Basis2D>
T
H1NormPreconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, long k1,
                                              XType xtype2, int j2, long k2) const
{
    T dd_x, id_x, dd_y, id_y;
    if (    (flens::IsSame<Basis_x, Basis<T,Primal,R,CDF> >::value)
             && (flens::IsSame<Basis_y, Basis<T,Primal,R,CDF> >::value) ){
        T tmp = 0.;
        if (xtype1==XBSpline) {
            dd_x = pow2i<T>(2*j1)*refval_dd_bspline;
            id_x = refval_id_bspline;
        }
        else {
            dd_x = pow2i<T>(2*j1)*refval_dd_wavelet;
            id_x = refval_id_wavelet;
        }
        if (xtype2==XBSpline) {
            dd_y = pow2i<T>(2*j2)*refval_dd_bspline;
            id_y = refval_id_bspline;
        }
        else {
            dd_y = pow2i<T>(2*j2)*refval_dd_wavelet;
            id_y = refval_id_wavelet;
        }
    }
    else {
        dd_x = _integral_x(j1,k1,xtype1,1,j1,k1,xtype1,1);
        id_x = _integral_x(j1,k1,xtype1,0,j1,k1,xtype1,0);
        dd_y = _integral_y(j2,k2,xtype2,1,j2,k2,xtype2,1);
        id_y = _integral_y(j2,k2,xtype2,0,j2,k2,xtype2,0);
    }
    return 1./std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y);
}

template <typename T, typename Basis2D>
T
H1NormPreconditioner2D<T,Basis2D>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k);
}

}   // namespace lawa

