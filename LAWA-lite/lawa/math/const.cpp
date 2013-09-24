#include <lawa/math/const.h>

namespace lawa {

double Const<double>::EQUALITY_EPS = 1e-15;
double Const<double>::SQRT2 = M_SQRT2;
double Const<double>::R_SQRT2 = 1./M_SQRT2;

//todo: revise!!
long double Const<long double>::EQUALITY_EPS = 1e-18;
long double Const<long double>::SQRT2 = M_SQRT2;
long double Const<long double>::R_SQRT2 = 1.L/M_SQRT2;

} // namespace lawa

