#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_TIMESTEPPINGERRORCALCULATION_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_TIMESTEPPINGERRORCALCULATION_H 1

namespace lawa {

/* Calculate Errors using trapezoidal rule for integrals 
 *
 */
 
// Error at time t in L_2(Omega)
template<typename T, typename Basis1D>
T
calculateL2Error(T t, const flens::DenseVector<flens::Array<T> >& u, T (*sol)(T,T), 
                const Basis1D& basis, const int J, const double deltaX = 1./128.);

// Error at time t in H^1(Omega)
template<typename T, typename Basis1D>
T
calculateH1Error(T t, const flens::DenseVector<flens::Array<T> >& u, T (*sol)(T,T), T (*dx_sol)(T,T),
                 const Basis1D& basis, const int J, const double deltaX = 1./128.);

// Error in L_2(0,T; L_2)
//      U(_, k) contains solution at timestep k
template<typename T, typename Basis1D>
T
calculateL2_L2_Error(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, T (*sol)(T,T), 
                     const Basis1D& basis, const int J, const double deltaT, const int K, 
                     const double deltaX = 1./128.);

// Error in L_2(0,T; H^1)
//      U(_, k) contains solution at timestep k
template<typename T, typename Basis1D>
T
calculateL2_H1_Error(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, 
                     T (*sol)(T,T), T (*dx_sol)(T,T), const Basis1D& basis, const int J, 
                     const double deltaT, const int K, const double deltaX = 1./128.);
    
} // namespace lawa

#include <lawa/methods/uniform/postprocessing/timesteppingerrorcalculation.tcc>


#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_TIMESTEPPINGERRORCALCULATION_H
