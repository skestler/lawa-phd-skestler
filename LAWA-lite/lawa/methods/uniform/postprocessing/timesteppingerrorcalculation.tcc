
namespace lawa {

template<typename T, typename Basis1D>
T
calculateL2Error(T t, const flens::DenseVector<flens::Array<T> >& u, T (*sol)(T,T), 
                const Basis1D& basis, const int J, const double deltaX)
{    
    T L2_error = 0.5*(evaluate(basis, J, u, 0., 0) - sol(t, 0.))*(evaluate(basis, J, u, 0., 0) - sol(t, 0.));
    for(double x = deltaX; x < 1.; x += deltaX){
        T u_approx = evaluate(basis, J, u, x, 0);
        L2_error += (u_approx -  sol(t, x)) * (u_approx - sol(t,x));
    }
    L2_error += 0.5*(evaluate(basis, J, u, 1., 0) - sol(t, 1.))*(evaluate(basis, J, u, 1., 0) - sol(t, 1.));
    
    L2_error *= deltaX;
    
    return sqrt(L2_error);
}

template<typename T, typename Basis1D>
T
calculateH1Error(T t, const flens::DenseVector<flens::Array<T> >& u, T (*sol)(T,T), T (*dx_sol)(T,T),
                 const Basis1D& basis, const int J, const double deltaX)
{    
    T L2_error = 0.5*(evaluate(basis, J, u, 0., 0) - sol(t, 0.))*(evaluate(basis, J, u, 0., 0) - sol(t, 0.));
    T dx_L2_error = 0.5*(evaluate(basis, J, u, 0., 1) - dx_sol(t, 0.))*(evaluate(basis, J, u, 0., 1) - dx_sol(t, 0.));
    for(double x = deltaX; x < 1.; x += deltaX){
        T u_approx = evaluate(basis, J, u, x, 0);
        T dx_u_approx = evaluate(basis, J, u, x, 1);
        L2_error += (u_approx -  sol(t, x)) * (u_approx - sol(t,x));
        dx_L2_error += (dx_u_approx -  dx_sol(t, x)) * (dx_u_approx - dx_sol(t,x));
    }
    L2_error += 0.5*(evaluate(basis, J, u, 1., 0) - sol(t, 1.))*(evaluate(basis, J, u, 1., 0) - sol(t, 1.));
    dx_L2_error += 0.5*(evaluate(basis, J, u, 1., 1) - dx_sol(t, 1.))*(evaluate(basis, J, u, 1., 1) - dx_sol(t, 1.));
    
    L2_error *= deltaX;
    dx_L2_error *= deltaX;
    
    return sqrt(L2_error + dx_L2_error);
}

template<typename T, typename Basis1D>
T
calculateL2_L2_Error(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, T (*sol)(T,T), 
                     const Basis1D& basis, const int J, const double deltaT, const int K, 
                     const double /*deltaX*/)
{
    T L2_error = 0;
    
    T factor, val_L2;
    for(int k = 0; k <= K; ++k){
        if ((k == 0) || (k == K)) {
            factor = 0.5;
        }
        else {
            factor = 1.;
        }
        flens::DenseVector<flens::Array<T> > u = U(U.rows(), k);
        val_L2 = calculateL2Error(k*deltaT, u, sol, basis, J, 1./pow2i<T>(J+2));
        L2_error += factor *  val_L2 * val_L2;
    }
    L2_error *= deltaT;
    
    return std::sqrt(L2_error);
}

template<typename T, typename Basis1D>
T
calculateL2_H1_Error(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, 
                     T (*sol)(T,T), T (*dx_sol)(T,T), const Basis1D& basis, const int J, 
                     const double deltaT, const int K, const double /*deltaX*/)
{
    T H1_error = 0;
    T factor, val_H1;
    for(int k = 0; k <= K; ++k){
        if ((k == 0) || (k == K)) {
            factor = 0.5;
        }
        else {
            factor = 1.;
        }
        flens::DenseVector<flens::Array<T> > u = U(U.rows(), k);
        val_H1 = calculateH1Error(k*deltaT, u, sol, dx_sol, basis, J, 1./pow2i<T>(J+2));
        H1_error += factor * val_H1 * val_H1;
    }
    H1_error *= deltaT;
    
    return std::sqrt(H1_error);
}
    
} // namespace lawa

