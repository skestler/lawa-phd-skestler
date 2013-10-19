namespace lawa {

template <typename T>
int
RefSols_PDE_RPlus1D<T>::nr;

template <typename T>
T
RefSols_PDE_RPlus1D<T>::reaction;

template <typename T>
T
RefSols_PDE_RPlus1D<T>::convection;

template <typename T>
T
RefSols_PDE_RPlus1D<T>::diffusion;

template <typename T>
flens::DenseVector<Array<T> >
RefSols_PDE_RPlus1D<T>::sing_pts;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
RefSols_PDE_RPlus1D<T>::deltas;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
RefSols_PDE_RPlus1D<T>::H1_deltas;


template <typename T>
void
RefSols_PDE_RPlus1D<T>::setExample(int _nr, T _reaction, T _convection, T _diffusion)
{
    nr=_nr;
    reaction      = _reaction;
    convection    = _convection;
    diffusion     = _diffusion;

    if (nr==1) {
        sing_pts.engine().resize(1);
        sing_pts(1) = 1./3.;
        deltas.engine().resize(1,2);
        deltas(1,1) = 1./3.; deltas(1,2) = diffusion*(1.1*exp(1./3.)-0.1);
        H1_deltas.engine().resize(1,2);
        H1_deltas(1,1) = 1./3.; H1_deltas(1,2) = 1.1*exp(1./3.)-0.1;
    }
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::exact(T x, int deriv)
{

    if      (nr==1) {
        if (deriv==0) {
            if (0<=x && x<1./3.) return (exp(x)-1.);
            else                 return exp(-0.1*(x-1./3.))*(exp(1./3.)-1.);
        }
        else if (deriv==1)  {
            if (0<=x && x<1./3.) return exp(x);
            else                 return -0.1*exp(-0.1*(x-1./3.))*(exp(1./3.)-1.);
        }
        else {
            if (0<=x && x<1./3.) return exp(x);
            else                 return 0.01*exp(-0.1*(x-1./3.))*(exp(1./3.)-1.);
        }
    }
    else {
        assert(0);
        return 0;
    }
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::u(T x)
{
    return RefSols_PDE_RPlus1D<T>::exact(x, 0);
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::d_u(T x)
{
    return RefSols_PDE_RPlus1D<T>::exact(x, 1);
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::rhs(T x)
{
    return -diffusion*exact(x,2) + convection*exact(x,1) + reaction*exact(x,0);
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::H1_rhs(T x)
{
    return -exact(x,2) + exact(x,0);
}

template <typename T>
T
RefSols_PDE_RPlus1D<T>::H1norm()
{
    if (nr==1) {
        long double L2norm_sq = (1./6.)*(41. - 72.*exp(1./3.) + 33.*exp(2./3.));
        long double H1semi_sq = (1./20.)*(-9.-2.*exp(1./3.) + 11.*exp(2./3.));
        return std::sqrt(L2norm_sq+H1semi_sq);
    }

    else {
        assert(0);
        return 0;
    }
}

} //namespace lawa
