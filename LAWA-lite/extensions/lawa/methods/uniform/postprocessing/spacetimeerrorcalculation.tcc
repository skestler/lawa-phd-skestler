#include <lawa/methods/uniform/datastructures/uniformindex2d.h>

namespace lawa{
    
template<typename T, typename Basis2D>
T
estimateSpaceTimeH1Error(const Basis2D& basis,
        const flens::DenseVector<flens::Array<T> >& u_approx,const int J_t_approx, const int J_x_approx, 
        const flens::DenseVector<flens::Array<T> >& u_exact, const int J_t_exact, const int J_x_exact)
{
    UniformIndex2D<Basis2D>  I_approx(basis, J_t_approx, J_x_approx);
    UniformIndex2D<Basis2D>  I_exact(basis, J_t_exact, J_x_exact);

    typename Basis2D::FirstBasisType b1 = basis.first;
    typename Basis2D::SecondBasisType b2 = basis.second;
    T error = 0.;

    // SF x SF
    for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
        for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
            int i_a = I_approx(XBSpline, b1.j0, kt, XBSpline, b2.j0, kx);
            int i_e = I_exact(XBSpline, b1.j0, kt, XBSpline, b2.j0, kx);
            error += (pow2i<T>(2*b1.j0 - 2*b2.j0) + pow2i<T>(2*b2.j0))
                     *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
        }
    }
    // SF x W
    int J_x_approx_max = basis.J2_max(J_t_approx, J_x_approx, b1.j0-1);
    int J_x_exact_max = basis.J2_max(J_t_exact, J_x_exact, b1.j0-1);
    for(int jx = b2.j0; jx <= std::min(J_x_approx_max, J_x_exact_max) - 1; ++jx){
        for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
            for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                int i_a = I_approx(XBSpline, b1.j0, kt, XWavelet, jx, kx);
                int i_e = I_exact(XBSpline, b1.j0, kt, XWavelet, jx, kx);
                error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                         *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
            }
        } 
    }
    for(int jx = std::min(J_x_approx_max, J_x_exact_max); jx <= std::max(J_x_approx_max, J_x_exact_max) - 1; ++jx){
        if(J_x_exact_max > J_x_approx_max){
            // approx coefficients are zero
            for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_e = I_exact(XBSpline, b1.j0, kt, XWavelet, jx, kx);
                    error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                             *u_exact(i_e)*u_exact(i_e);
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_a = I_approx(XBSpline, b1.j0, kt, XWavelet, jx, kx);
                    error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                             *u_approx(i_a)*u_approx(i_a);
                }
            }
        }
    }

    // W x SF
    int J_t_approx_max = basis.J1_max(J_t_approx, J_x_approx, b2.j0-1);
    int J_t_exact_max = basis.J1_max(J_t_exact, J_x_exact, b2.j0-1);
    for(int jt = b1.j0; jt <= std::min(J_t_approx_max, J_t_exact_max) - 1; ++jt){
        for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
            for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                int i_a = I_approx(XWavelet, jt, kt, XBSpline, b2.j0, kx);
                int i_e = I_exact(XWavelet, jt, kt, XBSpline, b2.j0, kx);
                error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                         *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
            }
        } 
    }
    for(int jt = std::min(J_t_approx_max, J_t_exact_max); jt <= std::max(J_t_approx_max, J_t_exact_max) - 1; ++jt){        
        if(J_t_exact_max > J_t_approx_max){
            // approx coefficients are zero
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                    int i_e = I_exact(XWavelet, jt, kt, XBSpline, b2.j0, kx);
                    error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                             *u_exact(i_e)*u_exact(i_e);
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                    int i_a = I_approx(XWavelet, jt, kt, XBSpline, b2.j0, kx);
                    error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                             *u_approx(i_a)*u_approx(i_a);
                }
            }
        }
    }

    // W x W    
    J_t_approx_max = basis.J1_max(J_t_approx, J_x_approx, b2.j0);
    J_t_exact_max = basis.J1_max(J_t_exact, J_x_exact, b2.j0);
    for(int jt = b1.j0; jt <= std::min(J_t_approx_max, J_t_exact_max) - 1; ++jt){
        J_x_approx_max = basis.J2_max(J_t_approx, J_x_approx, jt);
        J_x_exact_max = basis.J2_max(J_t_exact, J_x_exact, jt);
        for(int jx = b2.j0; jx <= std::min(J_x_approx_max, J_x_exact_max) - 1; ++jx){
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_a = I_approx(XWavelet, jt, kt, XWavelet, jx, kx);
                    int i_e = I_exact(XWavelet, jt, kt, XWavelet, jx, kx);
                    error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                             *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
                }
            }
        }
        for(int jx = std::min(J_x_approx_max, J_x_exact_max); jx <= std::max(J_x_approx_max, J_x_exact_max) - 1; ++jx){  
            if(J_t_exact_max > J_t_approx_max){
                // approx coefficients are zero
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_e = I_exact(XWavelet, jt, kt, XWavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_exact(i_e)*u_exact(i_e);
                    }
                }
            }
            else{
                // exact coefficients are zero
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_a = I_approx(XWavelet, jt, kt, XWavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_approx(i_a)*u_approx(i_a);
                    }
                }
            }
        } 
    }

    for(int jt = std::min(J_t_approx_max, J_t_exact_max); jt <= std::max(J_t_approx_max, J_t_exact_max) - 1; ++jt){          
        J_x_approx_max = basis.J2_max(J_t_approx, J_x_approx, jt);
        J_x_exact_max = basis.J2_max(J_t_exact, J_x_exact, jt);      
        if(J_t_exact_max > J_t_approx_max){
            // approx coefficients are zero
            for(int jx = b2.j0; jx <= J_x_exact_max - 1; ++jx){                                
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_e = I_exact(XWavelet, jt, kt, XWavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_exact(i_e)*u_exact(i_e);
                    }
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int jx = b2.j0; jx <= J_x_approx_max - 1; ++jx){                                
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_a = I_approx(XWavelet, jt, kt, XWavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_approx(i_a)*u_approx(i_a);
                    }
                }
            }
        }
    }


    return sqrt(error);
}

template<typename T, typename Basis2D>
T
calculateSpaceTimeL2Error(const Basis2D& basis, T (*sol)(T, T), T (*dx_sol)(T,T),
        const flens::DenseVector<flens::Array<T> > u, const int J_t, const int J_x, 
        const double deltaT, const double deltaX)
{
    T L2error = 0.;
    for(double t = 0.; t <= 1.; t += deltaT){
        T factor_t = deltaT;
        if((t == 0) || (t == 1.)){
            factor_t *= 0.5;
        }
        T space_H1error = 0;
        for(double x = 0; x <= 1.; x += deltaX){
            T factor_x = deltaX;
            if((x == 0) || (x == 1.)){
                factor_x *= 0.5;
            }
            T u_approx = evaluate(basis, J_t, J_x, u, t, x, 0, 0);
            T dx_u_approx = evaluate(basis, J_t, J_x, u, t, x, 0, 1);
            space_H1error += factor_x * ((u_approx -  sol(t, x)) * (u_approx - sol(t,x))
                                       + (dx_u_approx - dx_sol(t,x))*(dx_u_approx - dx_sol(t,x)));
        }
        L2error += factor_t * space_H1error;
    }
    
    return sqrt(L2error);
}                   
                       
                       
}                       

