namespace lawa {

template <typename T>
int
TensorRefSols_PDE_Interval_RPlus<T>::nr;

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::reaction;

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::convection_x;

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::convection_y;

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::diffusion_y;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Interval_RPlus<T>::sing_pts_x;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Interval_RPlus<T>::sing_pts_y;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Interval_RPlus<T>::deltas_x;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Interval_RPlus<T>::deltas_y;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Interval_RPlus<T>::H1_deltas_x;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Interval_RPlus<T>::H1_deltas_y;


template <typename T>
void
TensorRefSols_PDE_Interval_RPlus<T>::setExample(int _nr, T _reaction, T _convection_x, T _convection_y,
                                                T _diffusion_y)
{
    reaction     =_reaction;
    convection_x =_convection_x;
    convection_y =_convection_y;
    diffusion_y  =_diffusion_y;
    assert(reaction>=0);
    nr=_nr;

    if (nr==1) {
        sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
        sing_pts_y.engine().resize(1); sing_pts_y(1) = 1./3.;
        deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 1.5*exp(1./3.);
        deltas_y.engine().resize(1,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = diffusion_y*(7.*exp(1./3.)-6.);
        H1_deltas_x.engine().resize(1,2); H1_deltas_x(1,1) = 1./3.; H1_deltas_x(1,2) = 1.5*exp(1./3.);
        H1_deltas_y.engine().resize(1,2); H1_deltas_y(1,1) = 1./3.; H1_deltas_y(1,2) = (7.*exp(1./3.)-6.);;
    }
    else if (nr==2) {
        sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
        sing_pts_y.engine().resize(1); sing_pts_y(1) = 1./3.;
        deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 15*exp(1./3.);
        deltas_y.engine().resize(1,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = diffusion_y*(1.1*exp(1./3.)-0.1);
        H1_deltas_x.engine().resize(1,2); H1_deltas_x(1,1) = 1./3.; H1_deltas_x(1,2) = 15*exp(1./3.);
        H1_deltas_y.engine().resize(1,2); H1_deltas_y(1,1) = 1./3.; H1_deltas_y(1,2) = (1.1*exp(1./3.)-0.1);
    }
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact(T x, T y)
{
    return exact_x(x,0)*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_dx(T x, T y)
{
    return exact_x(x,1) * exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_dy(T x, T y)
{
    return exact_x(x,0) * exact_y(y,1);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_x(T x)
{
    return exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_y(T y)
{
    return exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::rhs_x(T x)
{
    return -exact_x(x,2) + convection_x*exact_x(x,1) + 0.5*reaction*exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::H1_rhs_x(T x)
{
    return -exact_x(x,2) + 0.5*exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::rhs_y(T y)
{
    return -diffusion_y*exact_y(y,2) + convection_y*exact_y(y,1) + 0.5*reaction*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::H1_rhs_y(T y)
{
    return -exact_y(y,2) + 0.5*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_x(T x, int deriv_x)
{
    if (nr==1) {
        if (deriv_x==0) {
            if (0<=x && x<1./3.) return exp(x)-1.;
            else                 return exp(-0.5*(x-1.))-1.;
        }
        else if (deriv_x==1)  {
            if (0<=x && x<1./3.) return exp(x);
            else                 return -0.5*exp(-0.5*(x-1.));
        }
        else {
            if (0<=x && x<1./3.) return exp(x);
            else                 return 0.25*exp(-0.5*(x-1.));
        }
    }
    else if (nr==2) {
        if (deriv_x==0) {
            if (0<=x && x<1./3.) return 10*exp(x)-10.;
            else                 return 10*exp(-0.5*(x-1.))-10.;
        }
        else if (deriv_x==1)  {
            if (0<=x && x<1./3.) return 10*exp(x);
            else                 return -5*exp(-0.5*(x-1.));
        }
        else {
            if (0<=x && x<1./3.) return 10*exp(x);
            else                 return 2.5*exp(-0.5*(x-1.));
        }
    }
    else {
        assert(0);
        return 0;
    }
    return 0;
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::exact_y(T y, int deriv_y)
{
    if (nr==1) {
        if (deriv_y==0) {
            if (0<=y && y<1./3.) return (exp(y)-1.);
            else                 return (1./std::pow(y+1.-1./3.,6.))*(exp(1./3.)-1);
        }
        else if (deriv_y==1)  {
            if (0<=y && y<1./3.) return exp(y);
            else                 return (-6./std::pow(y+1.-1./3.,7.))*(exp(1./3.)-1);
        }
        else {
            if (0<=y && y<1./3.) return exp(y);
            else                 return (42./std::pow(y+1.-1./3.,8.))*(exp(1./3.)-1);
        }
    }
    else if (nr==2) {
        if (deriv_y==0) {
            if (0<=y && y<1./3.) return (exp(y)-1.);
            else                 return exp(-0.1*(y-1./3.))*(exp(1./3.)-1.);
        }
        else if (deriv_y==1)  {
            if (0<=y && y<1./3.) return exp(y);
            else                 return -0.1*exp(-0.1*(y-1./3.))*(exp(1./3.)-1.);
        }
        else {
            if (0<=y && y<1./3.) return exp(y);
            else                 return 0.01*exp(-0.1*(y-1./3.))*(exp(1./3.)-1.);
        }
    }
    else {
        assert(0);
        return 0;
    }
    return 0;
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::H1norm()
{
    T ret = 0.;
    if (nr==1) {
        ret = 1.;
    }
    else if (nr==2) {
        long double L2norm_x_sq = 50.*(11.-12.*exp(1./3.)+3.*exp(2./3.));
        long double L2norm_y_sq = (1./6.)*(41. - 72.*exp(1./3.) + 33.*exp(2./3.));
        long double H1semi_x_sq = 75.*(-1.+exp(2./3.));
        long double H1semi_y_sq = (1./20.)*(-9.-2.*exp(1./3.) + 11.*exp(2./3.));
        long double H1norm = std::sqrt(L2norm_x_sq*L2norm_y_sq + L2norm_x_sq*H1semi_y_sq + H1semi_x_sq*L2norm_y_sq);
        ret = H1norm;
    }
    return ret;
}

template <typename T>
T
TensorRefSols_PDE_Interval_RPlus<T>::Energynorm()
{
    T ret = 0.;
    if (nr==2) {
        long double L2norm_x_sq = 50.*(11.-12.*exp(1./3.)+3.*exp(2./3.));
        long double L2norm_y_sq = (1./6.)*(41. - 72.*exp(1./3.) + 33.*exp(2./3.));
        long double H1semi_x_sq = 75.*(-1.+exp(2./3.));
        long double H1semi_y_sq = (1./20.)*(-9.-2.*exp(1./3.) + 11.*exp(2./3.));
        long double H1norm = std::sqrt(L2norm_x_sq*L2norm_y_sq + diffusion_y*L2norm_x_sq*H1semi_y_sq
                                       +reaction*H1semi_x_sq*L2norm_y_sq);
        ret = H1norm;
    }
    return ret;
}

}   // namespace lawa
