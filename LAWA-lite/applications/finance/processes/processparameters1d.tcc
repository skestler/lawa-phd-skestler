namespace lawa {

template <typename T>
ProcessParameters1D<T,CGMY>::ProcessParameters1D(T _r, T _k_C, T _k_G, T _k_M, T _k_Y)
: r(_r), k_C(_k_C), k_G(_k_G), k_M(_k_M), k_Y(_k_Y)
{

}

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,CGMY> &processparameters) {
    s << "_CGMY_r_" << processparameters.r << "_C_" << processparameters.k_C
      << "_G_" << processparameters.k_G << "_M_" << processparameters.k_M
      << "_Y_" << processparameters.k_Y;
    return s;
}

template <typename T>
ProcessParameters1D<T,CGMYe>::ProcessParameters1D(T _r, T _k_C, T _k_G, T _k_M, T _k_Y, T _sigma)
: r(_r), k_C(_k_C), k_G(_k_G), k_M(_k_M), k_Y(_k_Y), sigma(_sigma)
{

}

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,CGMYe> &processparameters) {
    s << "_CGMYe_r_" << processparameters.r << "_C_" << processparameters.k_C
      << "_G_" << processparameters.k_G << "_M_" << processparameters.k_M
      << "_Y_" << processparameters.k_Y << "_sigma_" << processparameters.sigma;
    return s;
}

template <typename T>
ProcessParameters1D<T,BlackScholes>::ProcessParameters1D(T _r, T _sigma)
: r(_r), sigma(_sigma)
{

}

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters1D<T,BlackScholes> &processparameters) {
    s << "_BS_r_" << processparameters.r << "_sigma_" << processparameters.sigma;
    return s;
}

}   // namespace lawa
