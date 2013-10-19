namespace lawa {

template <typename T>
int
TensorRefSols_PDE_Realline2D<T>::nr;

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::reaction;

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::convection_x;

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::convection_y;

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::diffusion_y;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Realline2D<T>::sing_pts_x;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Realline2D<T>::sing_pts_y;

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::singpt_ex2 = 1./3.;
//TensorRefSols_PDE_Realline2D<T>::singpt_ex2 = 0.0001;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline2D<T>::deltas_x;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline2D<T>::deltas_y;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline2D<T>::H1_deltas_x;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline2D<T>::H1_deltas_y;


template <typename T>
void
TensorRefSols_PDE_Realline2D<T>::setExample(int _nr, T _reaction, T _convection_x, T _convection_y,
                                            T _diffusion_y)
{
    reaction     =_reaction;
    convection_x =_convection_x;
    convection_y =_convection_y;
    diffusion_y =_diffusion_y;
    assert(reaction>=0);
    nr=_nr;

    std::cerr << "TensorRefSols_PDE_Realline2D<T>: nr = " << nr << ", reaction = " << reaction
              << ", convection_x = " << convection_x << ", convection_y = " << convection_y
              << ", diffusion_y = " << diffusion_y << std::endl;

    if (nr==2) {
        sing_pts_x.engine().resize(1); sing_pts_x(1) = singpt_ex2;
        sing_pts_y.engine().resize(1); sing_pts_y(1) = singpt_ex2;
        deltas_x.engine().resize(1,2); deltas_x(1,1) = singpt_ex2; deltas_x(1,2) = 0.2;
        deltas_y.engine().resize(1,2); deltas_y(1,1) = singpt_ex2; deltas_y(1,2) = diffusion_y*1.;
        H1_deltas_x.engine().resize(1,2); H1_deltas_x(1,1) = singpt_ex2; H1_deltas_x(1,2) = 0.2;
        H1_deltas_y.engine().resize(1,2); H1_deltas_y(1,1) = singpt_ex2; H1_deltas_y(1,2) = 1.;
    }
    if (nr==3) {
        sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
        deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 4.;
        H1_deltas_x.engine().resize(1,2); H1_deltas_x(1,1) = 1./3.; H1_deltas_x(1,2) = 4.;
    }
    if (nr==5) {
        sing_pts_x.engine().resize(3); sing_pts_x = 0., 1./3., 1.;
        deltas_x.engine().resize(3,2);
        deltas_x(1,1) = 0.;    deltas_x(1,2) = -4.;
        deltas_x(2,1) = 1./3.; deltas_x(2,2) =  4.*exp(1./3.)+2.*exp(-0.5*(1./3.-1.));
        deltas_x(3,1) = 1.;    deltas_x(3,2) = -2.;
        H1_deltas_x.engine().resize(3,2);
        H1_deltas_x(1,1) = 0.;    H1_deltas_x(1,2) = -4.;
        H1_deltas_x(2,1) = 1./3.; H1_deltas_x(2,2) = 4.*exp(1./3.)+2.*exp(-0.5*(1./3.-1.));
        H1_deltas_x(3,1) = 1.;    H1_deltas_x(3,2) = -2.;

        sing_pts_y.engine().resize(3); sing_pts_y = 0., 1./3., 1.;
        deltas_y.engine().resize(3,2);
        deltas_y(1,1) = 0.;    deltas_y(1,2) = diffusion_y*(-4.);
        deltas_y(2,1) = 1./3.; deltas_y(2,2) = diffusion_y*(4.*exp(1./3.)+2.*exp(-0.5*(1./3.-1.)));
        deltas_y(3,1) = 1.;    deltas_y(3,2) = diffusion_y*(-2.);
        H1_deltas_y.engine().resize(3,2);
        H1_deltas_y(1,1) = 0.;    H1_deltas_y(1,2) = -4.;
        H1_deltas_y(2,1) = 1./3.; H1_deltas_y(2,2) = 4.*exp(1./3.)+2.*exp(-0.5*(1./3.-1.));
        H1_deltas_y(3,1) = 1.;    H1_deltas_y(3,2) = -2.;

    }
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact(T x, T y)
{
    return exact_x(x,0)*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact_dx(T x, T y)
{
    return exact_x(x,1) * exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact_dy(T x, T y)
{
    return exact_x(x,0) * exact_y(y,1);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact_x(T x)
{
    return exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact_y(T y)
{
    return exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::rhs_x(T x)
{
    return -exact_x(x,2) + convection_x*exact_x(x,1) + 0.5*reaction*exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::H1_rhs_x(T x)
{
    return -exact_x(x,2) + 0.5*exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::rhs_y(T y)
{
    return -diffusion_y*exact_y(y,2) + convection_y*exact_y(y,1) + 0.5*reaction*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::H1_rhs_y(T y)
{
    return -exact_y(y,2) + 0.5*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::exact_x(T x, int deriv_x)
{
    if (nr==1) {
        if (deriv_x==0)         return  std::exp(-0.1*(x+0.1)*(x+0.1));
        else if (deriv_x==1)    return -std::exp(-0.1*(x+0.1)*(x+0.1)) * 2*0.1*(x+0.1);
        else                     return  std::exp(-0.1*(x+0.1)*(x+0.1))
                                        *(4*0.1*0.1*(x+0.1)*(x+0.1)-2*0.1);
    }
    else if (nr==2) {
        if (deriv_x==0)             return  std::exp(-0.1*fabs(x-singpt_ex2));
        else if (deriv_x==1) {
            if (x < singpt_ex2)     return   0.1*std::exp(-0.1*fabs(x-singpt_ex2));
            else                    return  -0.1*std::exp(-0.1*fabs(x-singpt_ex2));
        }
        else                        return   0.01*std::exp(-0.1*fabs(x-singpt_ex2));
    }
    else if (nr==3) {
        if (deriv_x==0)         return  std::exp(-2.*fabs(x-1./3.));
        else if (deriv_x==1) {
            if (x < 1./3.)      return   2.*std::exp(-2.*fabs(x-1./3.));
            else                return  -2.*std::exp(-2.*fabs(x-1./3.));
        }
        else {
            if (x < 1./3.)      return   4.*std::exp(-2.*fabs(x-1./3.));
            else                return   4.*std::exp(-2.*fabs(x-1./3.));
        }
    }
    else if (nr==5) {
        if (deriv_x==0) {
            if (x>0 && x<1./3.)        return   4.*(exp(x)-1);
            else if (x>1./3. && x<1.)  return   4.*(exp(-0.5*(x-1.))-1);
            else                       return   0.;
        }
        else if (deriv_x==1) {
            if (x>0 && x<1./3.)        return   4.*exp(x);
            else if (x>1./3. && x<1.)  return   -2.*exp(-0.5*(x-1.));
            else                       return   0.;
        }
        else {
            if (x>0 && x<1./3.)        return   4.*exp(x);
            else if (x>1./3. && x<1.)  return      exp(-0.5*(x-1.));
            else                       return   0.;
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
TensorRefSols_PDE_Realline2D<T>::exact_y(T y, int deriv_y)
{
    if (nr==1) {
        if (deriv_y==0)         return  std::exp(-0.5*(y-0.1)*(y-0.1));
        else if (deriv_y==1)    return -std::exp(-0.5*(y-0.1)*(y-0.1)) * 2*0.5*(y-0.1);
        else                    return  std::exp(-0.5*(y-0.1)*(y-0.1))
                                        * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
    }
    else if (nr==2) {
        if (deriv_y==0)         return  std::exp(-0.5*fabs(y-singpt_ex2));
        else if (deriv_y==1) {
            if (y < singpt_ex2) return   0.5*std::exp(-0.5*fabs(y-singpt_ex2));
            else                return  -0.5*std::exp(-0.5*fabs(y-singpt_ex2));
        }
        else                    return   0.25*std::exp(-0.5*fabs(y-singpt_ex2));
    }
    else if (nr==3) {
        if (deriv_y==0)         return std::exp(-0.1*(y-1./3.)*(y-1./3.));
        else if (deriv_y==1)    return -0.2*(y-1./3.)*std::exp(-0.1*(y-1./3.)*(y-1./3.));
        else                    return  (0.04*(y-1./3)*(y-1./3.)-0.2)
                                       *std::exp(-0.1*(y-1./3.)*(y-1./3.));
    }
    else if (nr==5) {
        if (deriv_y==0) {
            if (y>0 && y<1./3.)        return   4.*(exp(y)-1);
            else if (y>1./3. && y<1.)  return   4.*(exp(-0.5*(y-1.))-1);
            else                       return   0.;
        }
        else if (deriv_y==1) {
            if (y>0 && y<1./3.)        return   4.*exp(y);
            else if (y>1./3. && y<1.)  return   -2.*exp(-0.5*(y-1.));
            else                       return   0.;
        }
        else {
            if (y>0 && y<1./3.)        return   4.*exp(y);
            else if (y>1./3. && y<1.)  return      exp(-0.5*(y-1.));
            else                       return   0.;
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
TensorRefSols_PDE_Realline2D<T>::H1norm()
{
    T ret = 0.;
    if (nr==1) {
        //T a1 = -0.1;
        //T a2 = 0.1;
        T b1 = 0.1;
        T b2 = 0.5;
        T c1=1., c2= 1.;
        ret += c1*c1*std::sqrt(0.5*M_PI/(b1)) * c2*c2*std::sqrt(0.5*M_PI/(b2));
        ret += c1*c1*b1*b1*std::sqrt(0.5*M_PI)/std::pow(b1,(T)1.5)*c2*c2*std::sqrt(0.5*M_PI/(b2));
        ret += c1*c1*std::sqrt(0.5*M_PI/(b1))*c2*c2*b2*b2*std::sqrt(0.5*M_PI)/std::pow(b2,(T)1.5);
        ret = sqrt(ret);
    }
    if (nr==2) {
        T i11   = 10.;
        T i22   = 2.;
        T ddi11 = 0.1;
        T ddi22 = 0.5;
        ret = sqrt(ddi11*i22 + i11*ddi22 + i11*i22);
    }
    if (nr==3) {
        ret += 3.96332729760601*0.5 + 0.3963327297606012*0.5 + 3.96332729760601*2.;
        return sqrt(ret);
    }
    if (nr==5) {
        T i1 = 8.*(11.-12.*exp(1./3.)+3.*exp(2./3.));
        T i2 = 12.*(-1+exp(2./3.));
        return sqrt(i1*i1+2*i1*i2);
    }
    return ret;
}

template <typename T>
T
TensorRefSols_PDE_Realline2D<T>::Energynorm()
{
    T ret = 0.;
    if (nr==2) {
        T i11   = 10.;
        T i22   = 2.;
        T ddi11 = 0.1;
        T ddi22 = 0.5;
        ret = sqrt(ddi11*i22 + diffusion_y*i11*ddi22 + reaction*i11*i22);
    }
    return ret;
}

}   // namespace lawa
