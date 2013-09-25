namespace lawa {

template <typename T>
CharacteristicFunction1D<T,CGMY>::CharacteristicFunction1D
                                  (const ProcessParameters1D<T,CGMY> &_processparameters)
 : processparameters(_processparameters), drift(0.)
{
    T C = processparameters.k_C;
    T G = processparameters.k_G;
    T M = processparameters.k_M;
    T Y = processparameters.k_Y;
    drift = -C*boost::math::tgamma(-Y)*(std::pow(M-1,Y) - std::pow(M,Y) + std::pow(G+1,Y)
            - std::pow(G,Y));
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMY>::operator()(T t, T v)
{
    T C = processparameters.k_C;
    T G = processparameters.k_G;
    T M = processparameters.k_M;
    T Y = processparameters.k_Y;
    gsl_complex z1 = gsl_complex_rect(M,-v);
    gsl_complex z2 = gsl_complex_rect(G, v);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));
    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    tmp = gsl_complex_add(tmp, gsl_complex_rect(0.0,v*t*drift));
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMY>::operator()(T t, gsl_complex v)
{
    T C = processparameters.k_C;
    T G = processparameters.k_G;
    T M = processparameters.k_M;
    T Y = processparameters.k_Y;

    gsl_complex z1 = gsl_complex_rect(M+v.dat[1],-v.dat[0]);
    gsl_complex z2 = gsl_complex_rect(G-v.dat[1], v.dat[0]);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));

    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    tmp = gsl_complex_add(tmp, gsl_complex_rect(-v.dat[1]*t*drift,v.dat[0]*t*drift));
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMY>::zeta_0(void)
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
