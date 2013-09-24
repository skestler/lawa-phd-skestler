#ifndef LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H
#define LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H 1

#include <lawa/flensforlawa.h>

namespace lawa{
    
/* Fixed Point Solver for periodic problems
 *      This class solves periodic problems by repeatedly calling
 *      an underlying method (e.g. TimeStepping) with u_0^(i) = u(T)^(i-1)
 *      until || u_0^(i) - u(T)^(i) ||_l2 < tol 
 */
template<typename T, typename Method>
class FixedPointSolver
{
    public:
        typedef typename Method::RHSType RHSType;
        
        FixedPointSolver(Method& _method);
        
        flens::DenseVector<flens::Array<T> >
        solve(flens::DenseVector<flens::Array<T> > u_0, bool saveSols = false, 
              int maxIterations = 1000, T tol = 10e-15);
        
        flens::DenseVector<flens::Array<T> >
        solve(flens::DenseVector<flens::Array<T> > u_0,
                flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& fmatrix,
                int maxIterations = 1000, T tol = 10e-15);
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions(){ return method.getSolutions();}
    
        void setLevel(int level){ method.setLevel(level);}
        void setRHS(RHSType& rhs);
        
        
    private:
        Method& method;
        
        T
        getError(flens::DenseVector<flens::Array<T> >& u1, flens::DenseVector<flens::Array<T> >& u2);
    
};    
  
} // namespace lawa

#include <lawa/methods/uniform/solvers/fixedpointsolver.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H

