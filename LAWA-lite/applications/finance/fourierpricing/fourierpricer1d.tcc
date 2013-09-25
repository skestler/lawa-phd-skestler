namespace lawa {

template <typename T, ProcessType1D PType>
FourierPricer1D<T,PType>::FourierPricer1D(CharacteristicFunction1D<T,PType> &_charfunc,
                                          T _S0, T _maturity, T _K1, T _K2)
: charfunc(_charfunc), S0(_S0), maturity(_maturity), r(_charfunc.processparameters.r),
  K1(_K1), K2(_K2)
{
    assert(S0>0); assert(r>=0); assert(maturity>0); assert(K1>0); assert(K2>0);
}

/* This algorithm implements the fast Fourier transform for the computation of European call option
 * prices as described in [Cont & Tankov, 2004] (p.366). However, note that the return value is the
 * price of an European put option (obtained via put-call-parity).
 */
template <typename T, ProcessType1D PType>
void
FourierPricer1D<T,PType>::solve(T A, int J)
{
    int N = 1 << J;
    double delta = A/T(N-1);
    //T data[2*N];
    double *data;
    data = (double*)malloc( 2*N * sizeof(double) );
    double a=-std::log(K1/S0);

    for (int i = 0; i < N; i++)
    {
        double x = -0.5*A+delta*i;
        gsl_complex tmp = zeta(x);
        tmp = gsl_complex_mul(tmp,gsl_complex_exp(gsl_complex_rect(0,a*i*delta)));
        data[2*i] = tmp.dat[0]; data[2*i+1] = tmp.dat[1];
    }
    gsl_fft_complex_radix2_forward (data, 1, N);

    int n = std::ceil(((double)N/(double)(N-1)*A/(2*M_PI))*(a+log(K2/S0)))+1;
    //T _xa[n];
    //T _ya[n];
    double *_xa;
    double *_ya;
    _xa = (double*)malloc( n * sizeof(double) );
    _ya = (double*)malloc( n * sizeof(double) );


    for (int i = 0; i < n; i++)
    {
        gsl_complex fft = gsl_complex_rect(data[2*i],data[2*i+1]); //todo: tankov (p.366) funktionswerte abziehen
        double u = -a+(2*M_PI*i)/(N*delta);
        gsl_complex factor = gsl_complex_exp(gsl_complex_rect(0.0,u*0.5*A));
        gsl_complex res = gsl_complex_mul(fft,factor);

        double price = (A/(N*2.0*M_PI))*S0*res.dat[0]+S0*std::max(1.0-exp(u-r*maturity),0.0)-S0+exp(-r*maturity)*S0*exp(u);
        //double imag = (A/(N*2.0*M_PI))*res.dat[1];
        _xa[i] = S0*exp(u);
        _ya[i] = price;
        //std::cout << _xa[i] << " " << _ya[i] << std::endl;
    }
    interpol = gsl_interp_alloc(gsl_interp_linear,n);
    gsl_interp_init(interpol,_xa,_ya,n);
    accel = gsl_interp_accel_alloc();
    xa = _xa; ya =_ya;
    free(data);
    free(_xa);
    free(_ya);
    //std::cout << "FourierPricer: solve finished." << std::endl;
}

template <typename T, ProcessType1D PType>
T
FourierPricer1D<T,PType>::operator()(T K)
{
    return gsl_interp_eval(interpol,xa,ya,K,accel);
}

/* This algorithm implements the function zeta as defined in [Cont & Tankov, 2004](Eq. 11.17).
 */
template <typename T, ProcessType1D PType>
gsl_complex
FourierPricer1D<T,PType>::zeta(T v)
{
    gsl_complex res;
    if (fabs(v) > 10e-10) {
        gsl_complex z = gsl_complex_rect(v,-1);     // z = v-i
        gsl_complex charfunc_z = charfunc(maturity,z);
        gsl_complex zeta = gsl_complex_add_real(charfunc_z,-1.0);
        gsl_complex denom= gsl_complex_rect(-v*v,v);
        zeta = gsl_complex_div(zeta,denom);
        zeta = gsl_complex_mul(gsl_complex_exp(gsl_complex_rect(0,v*r*maturity)),zeta);
        res = zeta;
    }
    else  // v=0 -> 0/0 requires special treatment
    {
        res = gsl_complex_mul_real(charfunc.zeta_0(),maturity);
    }
    return res;
}

}   // namespace lawa
