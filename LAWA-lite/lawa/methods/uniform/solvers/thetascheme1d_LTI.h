#ifndef LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTI_H
#define LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTI_H 1

#include <lawa/settings/enum.h>
#include <lawa/integrals/integrals.h>
#include <lawa/operators/pdeoperators1d/pdeoperators1d.h>

namespace lawa{

/* ThetaScheme:
 *      This class solves an implicit linear system that arises in each time step of
 *      a time-stepping scheme for a linear and time-constant operator.
 *      It assumes a time-constant bilinear form, so that the system matrices
 *      are only assembled once. 
 *
 *      time_constant_rhs: rhsvector is only assembled once if true.
 *      use_pcg          : standard linear solver is gmres (false), otherwise cg (true)
 *      assemble_tol     : only save matrix entries with values > assemble_tol (absolute values)
 *      lintol           : solve linear system with accuracy lintol
 *      eta              : use weighted L2-scalar product \int w(x) v(x) e^{-2eta|x|} dx
 *      R1, R2           : use wavelets defined on [-R1, R2].
 *      order            : quadrature order for weighted L2-scalar product.
 */
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct=IdentityOperator1D<T,Basis> >
class ThetaScheme1D_LTI
{
    public: 
        typedef RHSIntegral RHSType;       
        
        ThetaScheme1D_LTI(const T _theta, const Basis& _basis, const BilinearForm& _a,
                          RHSIntegral& _rhs,
                          const bool time_constant_rhs=false, const bool _use_pcg=false,
                          T _assembletol=10e-15, T _lintol=10e-15);

        ThetaScheme1D_LTI(const T _theta, const Basis& _basis, const BilinearForm& _a,
                          RHSIntegral& _rhs, const L2ScalarProduct& _L2scalarproduct,
                          const bool time_constant_rhs=false, const bool _use_pcg=false,
                          T _assembletol=10e-15, T _lintol=10e-15);
    
        flens::DenseVector<flens::Array<T> > 
        solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, int level);
        
        flens::DenseVector<flens::Array<T> > 
        solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
              flens::DenseVector<flens::Array<T> > f, int level);
        
        void
        setRHS(RHSIntegral& _rhs);
        
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 
        getLHSMatrix(int level);                         
        
        // Adaptive Erweiterung: Timestep in jedem Lösungsschritt neu setzen,
        //flens::DenseVector<flens::Array<T> > 
        //solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level, T timestep);
        
        
    private:
        class Operator_LHSMatrix{
            private:
                ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>* scheme;
                const BilinearForm& a;
                T time_old;
                T time_new;
            
            public:                
                Operator_LHSMatrix(ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                                     L2ScalarProduct>* _scheme,
                                   const BilinearForm& _a);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;
                           
                void setTimes(T t1, T t2){ time_old = t1;
                                           time_new = t2;}
            
        };
        
        class Operator_RHSMatrix{
            private:
                const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                        L2ScalarProduct>* scheme;
                const BilinearForm& a;
                T time_old;
                T time_new;
            
            public:                
                Operator_RHSMatrix(const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                   L2ScalarProduct>* _scheme,
                                   const BilinearForm& _a);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;         
                       
                void setTimes(T t1, T t2){ time_old = t1;
                                           time_new = t2;}            
        };
        
        class Operator_RHSVector{
            private:
                const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                        L2ScalarProduct>* scheme;
                RHSIntegral& rhs;
                T time_old;
                T time_new;
                
            public:                
                Operator_RHSVector(const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                                           L2ScalarProduct>* _scheme,
                                   RHSIntegral& _rhs);
                
                T operator()(XType xtype, int j, int k) const;
                
                void setTimes(T t1, T t2){ time_old = t1;
                                           time_new = t2;}
                                           
                void setRHS(RHSIntegral& _rhs){ rhs = _rhs;}
                
            
        };
        
        friend class Operator_LHSMatrix;
        friend class Operator_RHSMatrix;
        friend class Operator_RHSVector;
        
        T theta;
        const Basis& basis;
        const IdentityOperator1D<T,Basis> standardL2scalarproduct;
        const L2ScalarProduct& L2scalarproduct;
        const bool time_constant_rhs;
        const bool use_pcg;
        T assembletol;
        T lintol;
        Assembler1D<T, Basis> assembler;

        //Integral<Gauss, Basis, Basis> integral;
        
        Operator_LHSMatrix op_LHSMatrix;
        Operator_RHSMatrix op_RHSMatrix;
        Operator_RHSVector op_RHSVector;
        
        DiagonalMatrixPreconditioner1D<T, Basis, Operator_LHSMatrix> prec;
        
        int currentLevel;
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > lhsmatrix;
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > rhsmatrix;
        flens::DiagonalMatrix<T>                                 P;
        flens::DenseVector<flens::Array<T> >                     rhsvector;
};
      
} // namespace lawa

#include <lawa/methods/uniform/solvers/thetascheme1d_LTI.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTI_H

