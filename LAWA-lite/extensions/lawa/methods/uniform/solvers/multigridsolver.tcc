namespace lawa {
    
    //typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver> 
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
MultigridSolver(PrimalBasis& _primalbasis, DualBasis& _dualbasis, Smoother& _smoother, 
                Solver& _solver, int _nu1, int _nu2, int _minLevel)
    : primalbasis(_primalbasis), dualbasis(_dualbasis), smoother(_smoother), solver(_solver),
      nu1(_nu1), nu2(_nu2), minLevel(_minLevel)
{
    //std::cout << "MultigridSolver: minLevel = " << minLevel << std::endl;
}

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
flens::DenseVector<flens::Array<T> >
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
vCycle(int i, int level, DenseVectorT& u, DenseVectorT& f)
{
    
    if(level == minLevel){
        solver.setLevel(level);
        u = solver.solve(u,f);
    }
    else{
        for(int j = 1; j <= i; ++j){
            // Smoothing
            smoother.setLevel(level);
            for(int i = 1; i <= nu1; ++i){
                u = smoother.solve(u,f);                
            }            
            // Residuum Restriction
            DenseVectorT r = smoother.getResiduum(u,f);
            DenseVectorT d;
            decompose(r, dualbasis, level-1, d);
            
            // Recursion
            DenseVectorT v(primalbasis.mra.cardI(level-1));
            DenseVectorT d_short = d(primalbasis.mra.rangeI(level-1));
            v = vCycle(1, level - 1, v, d_short);
            
            // Prolongation
            DenseVectorT v_long(dualbasis.mra_.rangeI_(level));
            v_long(primalbasis.mra.rangeI(level-1)) = v;
            DenseVectorT corr;
            reconstruct(v_long, primalbasis, level-1, corr);
            u = u + corr;
                        
            // Smoothing
            for(int i = 1; i <= nu2; ++i){
                u = smoother.solve(u,f);                
            }
        }
    }
    
    return u;
}

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
flens::DenseVector<flens::Array<T> >
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
wCycle(int i, int level, DenseVectorT& u, DenseVectorT& f)
{    
    if(level == minLevel){
        //std::cout << "Solver Level " << level << std::endl;
        solver.setLevel(level);
        u = solver.solve(u,f);
    }
    else{
        for(int j = 1; j <= i; ++j){
            // Smoothing
            //std::cout << "Smoothing Level " << level << std::endl;
            smoother.setLevel(level);
            for(int i = 1; i <= nu1; ++i){
                u = smoother.solve(u,f);                
            }            
            // Residuum Restriction
            //std::cout << "Restriction Level " << level << std::endl;
            DenseVectorT r = smoother.getResiduum(u,f);
            DenseVectorT d;
            decompose(r, dualbasis, level-1, d);
            
            // Recursion
            DenseVectorT v(primalbasis.mra.cardI(level-1));
            DenseVectorT d_short = d(primalbasis.mra.rangeI(level-1));
            //std::cout << "Recursion Level " << level << std::endl;
            v = wCycle(2, level - 1, v, d_short);
            
            // Prolongation
            DenseVectorT v_long(dualbasis.mra_.rangeI_(level));
            v_long(primalbasis.mra.rangeI(level-1)) = v;
            DenseVectorT corr;
            //std::cout << "Prolongation Level " << level << std::endl;
            reconstruct(v_long, primalbasis, level-1, corr);
            u = u + corr;
            
            
            // Smoothing
            for(int i = 1; i <= nu2; ++i){
                u = smoother.solve(u,f);                
            }
        }
    }
    
    return u;
}
}

