#ifndef LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H
#define LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H 1

#include <extensions/extensions.h>

namespace lawa{    

template<typename T, typename Basis>
class Assembler1D
{
    private:
        const Basis& basis;
   
    public: 
        Assembler1D(const Basis& _basis);

        /* Assemble Stiffness Matrix in transposed form, i.e.
         * corresponding to a(v,u)
         */
        template <typename BilinearForm>
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
        assembleStiffnessMatrix(BilinearForm& a, int J, T tol = 10e-15);
    
        template <typename RHSIntegral>
        flens::DenseVector<flens::Array<T> >
        assembleRHS(RHSIntegral& rhs, int J);
        
        template <typename Preconditioner>
        flens::DiagonalMatrix<T>    
        assemblePreconditioner(Preconditioner& P, int J);
};

} // namespace lawa

#include <lawa/methods/uniform/algorithms/assembler1d.tcc>

#endif // LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H

