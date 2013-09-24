namespace lawa{

// THETASCHEME
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
ThetaScheme1D_LTI(const T _theta, const Basis& _basis, const BilinearForm& _a, RHSIntegral& _rhs,
                  const bool _time_constant_rhs,
                  const bool _use_pcg, T _assembletol, T _lintol)
    : theta(_theta), basis(_basis),
      standardL2scalarproduct(basis), L2scalarproduct(standardL2scalarproduct),
      time_constant_rhs(_time_constant_rhs), use_pcg(_use_pcg),
      assembletol(_assembletol), lintol(_lintol),
      assembler(basis),
      //integral(_basis, _basis),
      op_LHSMatrix(this, _a), op_RHSMatrix(this, _a), op_RHSVector(this, _rhs), prec(op_LHSMatrix),
      currentLevel(-1), P(assembler.assemblePreconditioner(prec, basis.j0))
{
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
ThetaScheme1D_LTI(const T _theta, const Basis& _basis, const BilinearForm& _a, RHSIntegral& _rhs,
                  const L2ScalarProduct& _L2scalarproduct, const bool _time_constant_rhs,
                  const bool _use_pcg, T _assembletol, T _lintol)
    : theta(_theta), basis(_basis),
      standardL2scalarproduct(basis), L2scalarproduct(_L2scalarproduct),
      time_constant_rhs(_time_constant_rhs), use_pcg(_use_pcg),
      assembletol(_assembletol), lintol(_lintol),
      assembler(basis),
      //integral(_basis, _basis),
      op_LHSMatrix(this, _a), op_RHSMatrix(this, _a), op_RHSVector(this, _rhs), prec(op_LHSMatrix),
      currentLevel(-1), P(assembler.assemblePreconditioner(prec, basis.j0))
{
}
 
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
flens::DenseVector<flens::Array<T> > 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, int level)
{   
    op_RHSVector.setTimes(time_old, time_new);
    if(level != currentLevel){
        op_LHSMatrix.setTimes(time_old, time_new);
        op_RHSMatrix.setTimes(time_old, time_new);
        lhsmatrix = assembler.assembleStiffnessMatrix(op_LHSMatrix, level, assembletol);
        rhsmatrix = assembler.assembleStiffnessMatrix(op_RHSMatrix, level, assembletol);
        P = assembler.assemblePreconditioner(prec, level);
        rhsvector = assembler.assembleRHS(op_RHSVector, level);
        currentLevel = level;
    }
    if (!time_constant_rhs) {
        rhsvector = assembler.assembleRHS(op_RHSVector, level);
    }
    //std::cerr << "u_init    = " << u_init << std::endl;
    //std::cerr << "rhsvector = " << rhsvector << std::endl;
    flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + rhsvector;
    //std::cerr << "rhs       = " << rhs << std::endl;
    flens::DenseVector<flens::Array<T> > u(basis.mra.rangeI(level));
    //std::cerr << "First u   = " << u << std::endl;
    if (use_pcg) {
        pcg(P,lhsmatrix, u, rhs, lintol);
    }
    else {
        //int iters = gmres(lhsmatrix, u, rhs, lintol);
        int iters = pgmresm(P,lhsmatrix, u, rhs, lintol,10);
        //std::cerr << "Solve by gmres with iters = " << std::endl;
        //pgmresm(P,lhsmatrix, u, rhs, lintol, 20);
    }
    //std::cerr << "Next u    = " << u << std::endl;
    //std::cout << cg(lhsmatrix, u, rhs) << "cg iterations" << std::endl;
    //std::cout << pcg(P, lhsmatrix, u, rhs) << "pcg iterations" << std::endl;
    
    //std::cout << "u(" << time_new << "): " << u << std::endl; 
    return u;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
flens::DenseVector<flens::Array<T> > 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
      flens::DenseVector<flens::Array<T> > f, int level)
{
     if(level != currentLevel){
         op_LHSMatrix.setTimes(time_old, time_new);
         op_RHSMatrix.setTimes(time_old, time_new);
         lhsmatrix = assembler.assembleStiffnessMatrix(op_LHSMatrix, level, assembletol);
         rhsmatrix = assembler.assembleStiffnessMatrix(op_RHSMatrix, level, assembletol);
         P = assembler.assemblePreconditioner(prec, level);
         currentLevel = level;
     }
     flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + f;
     flens::DenseVector<flens::Array<T> > u(u_init);
     if (use_pcg) {
         pcg(P,lhsmatrix, u, rhs, lintol);
     }
     else {
         pgmres(P,lhsmatrix, u, rhs, lintol);
     }
     //std::cout << cg(lhsmatrix, u, rhs) << "cg iterations" << std::endl;
     //std::cout << pcg(P, lhsmatrix, u, rhs) << "pcg iterations" << std::endl;
     //std::cout << "u(" << time_new << "): " << u << std::endl; 
     return u;     
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
void
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
setRHS(RHSIntegral& _rhs)
{
    op_RHSVector.setRHS(_rhs);
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::
getLHSMatrix(int level)
{   
    if (level != currentLevel) {
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > matrix
                    = assembler.assembleStiffnessMatrix(op_LHSMatrix, level);
        return matrix;
    }
    return lhsmatrix;
}

/*======================================================================================*/    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,L2ScalarProduct>::Operator_LHSMatrix::
Operator_LHSMatrix(ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>* _scheme,
                   const BilinearForm& _a)
    : a(_a)
{   
    scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::Operator_LHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
    // (M + deltaT * theta * A_k+1)
    //return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0)
    //        + (time_new - time_old) * scheme->theta * a(xtype1, j1, k1, xtype2, j2, k2);
    return scheme->L2scalarproduct(xtype1, j1, k1, xtype2, j2, k2)
              + (time_new - time_old) * scheme->theta * a(xtype1, j1, k1, xtype2, j2, k2);
}


/*======================================================================================*/    
// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSMatrix::
Operator_RHSMatrix(const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                           L2ScalarProduct>* _scheme,
                   const BilinearForm& _a)
    : a(_a)
{
     scheme = _scheme;
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
   // (M - deltaT * (1-theta) * A_k)
   //return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0)
   //     - (time_new - time_old) * (1. - scheme->theta) * a(xtype1,j1,k1, xtype2, j2,k2);
    return scheme->L2scalarproduct(xtype1, j1, k1, xtype2, j2, k2)
          - (time_new - time_old) * (1. - scheme->theta) * a(xtype1,j1,k1, xtype2, j2,k2);
}

/*======================================================================================*/    
// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSVector::
Operator_RHSVector(const ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral,
                                           L2ScalarProduct>* _scheme,
                   RHSIntegral& _rhs)
    : rhs(_rhs)
{
     scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI<T, Basis, BilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSVector::
operator()(XType xtype, int j, int k) const
{   
    // deltaT * (theta*f_k+1 - (1-theta)*f_k)
    return (time_new - time_old)*(scheme->theta * rhs(time_new, xtype, j, k)
                    + (1. - scheme->theta)*rhs(time_old, xtype, j, k));
} 
  
}

