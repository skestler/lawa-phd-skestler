#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H 1

namespace lawa {

template<typename T, typename Basis1D>
void 
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis1D& basis, const int J, 
        const char* filename, const double deltaX=1./128.);
        
template<typename T, typename Basis2D>
void
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis2D& basis, const int J_x, const int J_y, 
                   const char* filename, const double deltaX=0.01, const double deltaY=0.01);
        
template<typename T, typename Basis1D>
void 
print_U(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, const Basis1D& basis, 
        const int J, const char* filename, const T timestep, const int K, const double deltaX=1./128.);
    
} // namespace lawa

#include <lawa/methods/uniform/postprocessing/plotting.tcc>

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H
