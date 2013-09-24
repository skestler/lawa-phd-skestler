#ifndef LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H
#define LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H 1

namespace lawa{

/* Multigrid Solver
 *      This class provides the framework for multigrid methods by
 *      implementing a general V- and W-Cycle for Wavelets 
 *      (using decompose/reconstruct for projections and restrictions) 
 */
template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
class MultigridSolver{
    
      typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
    
  private:
      PrimalBasis& primalbasis;
      DualBasis& dualbasis;
      Smoother& smoother;
      Solver& solver; 
      int nu1, nu2;                 
      int minLevel;
    
  public:
      MultigridSolver(PrimalBasis& _primalbasis, DualBasis& _dualbasis, Smoother& _smoother, 
                      Solver& _solver, int _nu1, int _nu2, int _minLevel = 0);
      
      DenseVectorT
      vCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);
                  
      DenseVectorT
      wCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);      
      
      int getMinLevel(){ return minLevel;}
};

} //  namespace lawa

#include <lawa/methods/uniform/solvers/multigridsolver.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H

