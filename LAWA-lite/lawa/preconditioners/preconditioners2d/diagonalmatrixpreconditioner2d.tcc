#include <cmath>

namespace lawa {

template <typename T, typename Basis2D, typename BilinearForm>
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::DiagonalMatrixPreconditioner2D(BilinearForm &a)
    : _a(a), theta(0.), timestep(0.)
{
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(XType xtype1, int j1, long k1,
                                                                   XType xtype2, int j2, long k2) const
{
    T val = fabs(_a(xtype1,j1,k1,xtype2,j2,k2,xtype1,j1,k1,xtype2,j2,k2));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(const Index2D &index) const
{
    T val = fabs(_a(index,index));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator[](const Index2D &index)
{
    T val = fabs(_a(index,index));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename Basis2D, typename BilinearForm>
void
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>
::setThetaTimeStepParameters(T _theta, T _timestep)
{
    theta = _theta;
    timestep = _timestep;
}

}   // namespace lawa

