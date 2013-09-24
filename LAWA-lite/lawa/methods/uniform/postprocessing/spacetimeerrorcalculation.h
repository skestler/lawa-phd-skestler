#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEERRORCALCULATION_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEERRORCALCULATION_H 1

namespace lawa{

/* estimateSpaceTimeH1error:
 * estimate error in H1(0,T; V') ^ L2(0,T; V)
 * using norm equivalence relations and the wavelet
 * representation of an "exact" solution.
 */
template<typename T, typename Basis2D>
T
estimateSpaceTimeH1Error(const Basis2D& basis,
        const flens::DenseVector<flens::Array<T> >& u_approx,const int J_t_approx, const int J_x_approx, 
        const flens::DenseVector<flens::Array<T> >& u_exact, const int J_t_exact, const int J_x_exact);

template<typename T, typename Basis2D>
T
calculateSpaceTimeL2Error(const Basis2D& basis, T (*sol)(T, T), T (*dx_sol)(T,T),
        const flens::DenseVector<flens::Array<T> > u, const int J_t, const int J_x, 
        const double deltaT=1./128., const double deltaX=1./128.);


} // namespace lawa

#include <lawa/methods/uniform/postprocessing/spacetimeerrorcalculation.tcc>

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEERRORCALCULATION_H

