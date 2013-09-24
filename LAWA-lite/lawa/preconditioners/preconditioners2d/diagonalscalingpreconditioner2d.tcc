#include <cassert>

namespace lawa {

template<typename T>
DiagonalScalingPreconditioner2D<T>::DiagonalScalingPreconditioner2D(int sx, int sy)
    : _sx(sx), _sy(sy)
{
    assert(sx >= 0);
    assert(sy >= 0);
}

template<typename T>
T
DiagonalScalingPreconditioner2D<T>::operator()(XType XisSpline, int jx, long /*kx*/,
                                               XType YisSpline, int jy, long /*ky*/) const
{
    assert(XisSpline);
    assert(YisSpline);

    return pow2i<T>(- _sx*jx - _sy*jy);
}

} // namespace lawa

