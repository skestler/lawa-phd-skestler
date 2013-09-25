namespace lawa {

template <typename T>
CharacteristicFunction1D<T,CGMYe>::CharacteristicFunction1D
                                   (const ProcessParameters1D<T,CGMYe> &_processparameters)
 : processparameters(_processparameters), drift(0.)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;
    drift = -C*boost::math::tgamma(-Y)*(std::pow(M-1,Y) - std::pow(M,Y) + std::pow(G+1,Y)
            - std::pow(G,Y))  - 0.5*sigma*sigma;;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::operator()(T t, T v)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;

    gsl_complex z1 = gsl_complex_rect(M,-v);
    gsl_complex z2 = gsl_complex_rect(G, v);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));
    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    tmp = gsl_complex_add(tmp, gsl_complex_rect(0.0,v*t*drift));
    tmp = gsl_complex_add_real(tmp,-0.5*sigma*sigma*v*v);
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::operator()(T t, gsl_complex v)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;

    gsl_complex z1 = gsl_complex_rect(M+v.dat[1],-v.dat[0]);
    gsl_complex z2 = gsl_complex_rect(G-v.dat[1], v.dat[0]);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    //gsl_complex z3 = gsl_complex_rect(-0.5*sigma*sigma*v.dat[0]*v.dat[0],
    //                                  +0.5*sigma*sigma*v.dat[1]*v.dat[1]);


    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    //tmp = gsl_complex_add(tmp,z3);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));

    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    gsl_complex z3 = gsl_complex_rect(t*0.5*sigma*sigma*(-v.dat[0]*v.dat[0]+v.dat[1]*v.dat[1]),
                                           -t*sigma*sigma*v.dat[0]*v.dat[1]);

    tmp = gsl_complex_add(tmp, gsl_complex_rect(-v.dat[1]*t*drift,v.dat[0]*t*drift));
    tmp = gsl_complex_add(tmp, z3);
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::zeta_0(void)
{
    T C = processparameters.k_C;
    T G = processparameters.k_G;
    T M = processparameters.k_M;
    T Y = processparameters.k_Y;
    gsl_complex res = gsl_complex_rect(C*boost::math::tgamma(-Y)*(-Y*std::pow(M-1,Y-1)
                                       +Y*std::pow(G+1,Y-1))+drift,0.);
    return res;
}

}   // namespace lawa
