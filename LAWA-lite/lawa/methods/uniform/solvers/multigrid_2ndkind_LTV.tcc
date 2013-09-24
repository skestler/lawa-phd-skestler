namespace lawa{

template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::
MultiGrid_2ndKind_LTV(PrimalBasis& _b, DualBasis& _b_, BilinearForm& a, RHSIntegral& rhs, T theta,
                  T deltaT, int timesteps, int minLevel)
    : b(_b), b_(_b_), full_theta(theta, b, a, rhs), full_ts(full_theta, deltaT, timesteps, 0),
      hom_theta(theta, b, a, hom_rhs), hom_ts(hom_theta, deltaT, timesteps, 0),
      fp(hom_ts), mg_smoother(this), mg_solver(this), mg(b, b_, mg_smoother, mg_solver, 1, 0, minLevel)
{ 
}

template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::
solve(flens::DenseVector<flens::Array<T> >& u0, int maxLevel)
{
    flens::DenseVector<flens::Array<T> > u(u0), zeros;
    for(int i = mg.getMinLevel() + 1; i <= maxLevel; ++i){
        flens::DenseVector<flens::Array<T> > u_prol;
        flens::DenseVector<flens::Array<T> > u_long(b_.mra_.rangeI_(i));
        u_long(b.mra.rangeI(i-1)) = u;
        reconstruct(u_long, b, i-1, u_prol);
        zeros.engine().resize(b.mra.cardI(i), b.mra.rangeI(i).firstIndex());
        full_ts.setLevel(i);
        flens::DenseVector<flens::Array<T> > f = full_ts.solve(zeros);                  
        u = mg.wCycle(1, i, u_prol, f);
    }
    
    return u;
}



template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::MG_2ndKind_Smoother::
MG_2ndKind_Smoother(MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref)
{
    mg_ptr = ref;
}

template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> >  
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::MG_2ndKind_Smoother::
solve(flens::DenseVector<flens::Array<T> > u, flens::DenseVector<flens::Array<T> > f)
{
   flens::DenseVector<flens::Array<T> > sol = mg_ptr->hom_ts.solve(u);
   return sol + f;
}



template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::MG_2ndKind_Solver::
MG_2ndKind_Solver(MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref)
{
    mg_ptr = ref;
}

template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
MultiGrid_2ndKind_LTV<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>::MG_2ndKind_Solver::
solve(flens::DenseVector<flens::Array<T> > u0, flens::DenseVector<flens::Array<T> > f)
{
    int steps = mg_ptr->hom_ts.getSteps();
    int deltaT = mg_ptr->hom_ts.getDeltaT();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > fmatrix(u0.length(), steps);
    
    
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
    lhsmatrix = mg_ptr->hom_theta.getLHSMatrix((steps-1)*deltaT, steps*deltaT, mg_ptr->mg.getMinLevel());
    
    flens::DenseVector<flens::Array<T> > rhsvector = lhsmatrix * f;
    fmatrix(flens::_, steps) = rhsvector;
    
    flens::DenseVector<flens::Array<T> > u = mg_ptr->fp.solve(u0, fmatrix);
    
    return u;
}

}

