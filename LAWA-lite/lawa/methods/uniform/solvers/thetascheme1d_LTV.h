#ifndef LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTV_H
#define LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTV_H 1

#include <lawa/settings/enum.h>
#include <lawa/preconditioners/preconditioners1d/diagonalmatrixpreconditioner1d.h>

namespace lawa{
    
/* ThetaScheme:
 *      This class solves an implicit linear system that arises in each time step of a 
 *      timestepping scheme for a linear but time-dependent operator.
 */
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
class ThetaScheme1D_LTV
{
    public: 
        typedef RHSIntegral RHSType;       
        
        ThetaScheme1D_LTV(const T _theta, const Basis& _basis, const BilinearForm& _a, RHSIntegral& _rhs);
    
        flens::DenseVector<flens::Array<T> > 
        solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, int level);
        
        flens::DenseVector<flens::Array<T> > 
        solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
              flens::DenseVector<flens::Array<T> > f, int level);
        
        void
        setRHS(RHSIntegral& _rhs);
        
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 
        getLHSMatrix(T time_old, T time_new, int level);
                           
        // Adaptive Erweiterung: Timestep in jedem Lösungsschritt neu setzen,
        //flens::DenseVector<flens::Array<T> > 
        //solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level, T timestep);
        
        
    private:
        class Operator_LHSMatrix{
            private:
                ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* scheme;
                const BilinearForm& a;
                T time_old;
                T time_new;
            
            public:                
                Operator_LHSMatrix(ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
                                   const BilinearForm& _a);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;

                
                void setTimes(T t1, T t2){ time_old = t1;
                                           time_new = t2;}
            
        };
        
        class Operator_RHSMatrix{
            private:
                const ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* scheme; 
                const BilinearForm& a;
                T time_old;
                T time_new;
            
            public:                
                Operator_RHSMatrix(const ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
                                   const BilinearForm& _a);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;
                
                void setTimes(T t1, T t2){ time_old = t1;
                                           time_new = t2;}            
            
        };
        
        class Operator_RHSVector{
            private:
                const ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* scheme; 
                RHSIntegral& rhs;
                T time_old;
                T time_new;
                
            public:                
                Operator_RHSVector(const ThetaScheme1D_LTV<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
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
        Assembler1D<T, Basis> assembler;
        
        Integral<Gauss, Basis, Basis> integral;
        
        Operator_LHSMatrix op_LHSMatrix;
        Operator_RHSMatrix op_RHSMatrix;
        Operator_RHSVector op_RHSVector;
        
        DiagonalMatrixPreconditioner1D<T, Basis, Operator_LHSMatrix> prec;
};
      
} // namespace lawa

#include <lawa/methods/uniform/solvers/thetascheme1d_LTV.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_THETASCHEME1D_LTV_H

