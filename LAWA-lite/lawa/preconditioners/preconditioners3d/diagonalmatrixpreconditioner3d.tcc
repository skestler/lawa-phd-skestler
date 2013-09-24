#include <cmath>

namespace lawa {

template <typename T, typename Basis3D, typename BilinearForm>
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::DiagonalMatrixPreconditioner3D(const BilinearForm &a)
    : _a(a)
{
}

template <typename T, typename Basis3D, typename BilinearForm>
T
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::operator()(XType xtype1, int j1, long k1,
                                                                   XType xtype2, int j2, long k2,
                                                                   XType xtype3, int j3, long k3) const
{
    return 1./std::sqrt(fabs(_a(xtype1,j1,k1, xtype2,j2,k2, xtype3,j3,k3,
                                xtype1,j1,k1, xtype2,j2,k2, xtype3,j3,k3)));
}

template <typename T, typename Basis3D, typename BilinearForm>
T
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::operator()(const Index3D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k,
                            index.index3.xtype, index.index3.j, index.index3.k);
}

} // namespace lawa

