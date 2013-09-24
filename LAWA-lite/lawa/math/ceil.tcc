#include <cmath>

namespace lawa {

template <typename T>
T
iceil(double x)
{
    return static_cast<T>(std::ceil(x));
}

} // namespace lawa
