namespace lawa {

template <typename T>
int
TensorRefSols_PDE_Realline3D<T>::nr;

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::c;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Realline3D<T>::sing_pts_x;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Realline3D<T>::sing_pts_y;

template <typename T>
flens::DenseVector<Array<T> >
TensorRefSols_PDE_Realline3D<T>::sing_pts_z;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline3D<T>::deltas_x;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline3D<T>::deltas_y;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorRefSols_PDE_Realline3D<T>::deltas_z;


template <typename T>
void
TensorRefSols_PDE_Realline3D<T>::setExample(int _nr, T _c)
{
    c=_c;
    assert(c>=0);
    nr=_nr;

    if (nr==2) {
        sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
        deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 4.;
        sing_pts_y.engine().resize(1); sing_pts_y(1) = 1./3.;
        deltas_y.engine().resize(1,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = 1.;
        sing_pts_z.engine().resize(1); sing_pts_z(1) = 1./3.;
        deltas_z.engine().resize(1,2); deltas_z(1,1) = 1./3.; deltas_z(1,2) = 0.2;
    }
    else if (nr==3) {
        sing_pts_z.engine().resize(1); sing_pts_z(1) = 1./3.;
        deltas_z.engine().resize(1,2); deltas_z(1,1) = 1./3.; deltas_z(1,2) = 4.;
    }
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact(T x, T y, T z)
{
    return exact_x(x,0)*exact_y(y,0)*exact_z(z,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_dx(T x, T y, T z)
{
    return exact_x(x,1) * exact_y(y,0) * exact_z(z,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_dy(T x, T y, T z)
{
    return exact_x(x,0) * exact_y(y,1) * exact_z(z,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_dz(T x, T y, T z)
{
    return exact_x(x,0) * exact_y(y,0) * exact_z(z,1);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_x(T x)
{
    return exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_y(T y)
{
    return exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_z(T z)
{
    return exact_z(z,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::rhs_x(T x)
{
    return -exact_x(x,2) + (1./3.)*c*exact_x(x,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::rhs_y(T y)
{
    return -exact_y(y,2) + (1./3.)*c*exact_y(y,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::rhs_z(T z)
{
    return -exact_z(z,2) + (1./3.)*c*exact_z(z,0);
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_x(T x, int deriv_x)
{
    if (nr==1) {
        if (deriv_x==0)         return  std::exp(-0.1*(x+0.1)*(x+0.1));
        else if (deriv_x==1)    return -std::exp(-0.1*(x+0.1)*(x+0.1)) * 2*0.1*(x+0.1);
        else                    return  std::exp(-0.1*(x+0.1)*(x+0.1))
                                        * (4*0.1*0.1*(x+0.1)*(x+0.1)-2*0.1);
    }
    else if (nr==2) {
        if (deriv_x==0)         return  std::exp(-2.*fabs(x-1./3.));
        else if (deriv_x==1) {
            if (x < 1./3.)      return   2.*std::exp(-2.*fabs(x-1./3.));
            else                return -2.*std::exp(-2.*fabs(x-1./3.));
        }
        else                    return   4.*std::exp(-2.*fabs(x-1./3.));
    }
    else if (nr==3) {
        if (deriv_x==0)         return  std::exp(-0.1*(x+0.1)*(x+0.1));
        else if (deriv_x==1)    return -std::exp(-0.1*(x+0.1)*(x+0.1)) * 2*0.1*(x+0.1);
        else                    return std::exp(-0.1*(x+0.1)*(x+0.1))
                                        * (4*0.1*0.1*(x+0.1)*(x+0.1)-2*0.1);
    }
    else {
        assert(0); return 0;
    }
    return 0;
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_y(T y, int deriv_y)
{
    if (nr==1) {
        if (deriv_y==0)         return  std::exp(-0.5*(y-0.1)*(y-0.1));
        else if (deriv_y==1)    return -std::exp(-0.5*(y-0.1)*(y-0.1)) * 2*0.5*(y-0.1);
        else                    return  std::exp(-0.5*(y-0.1)*(y-0.1))
                                        * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
    }
    else if (nr==2) {
        if (deriv_y==0)         return  std::exp(-0.5*fabs(y-1./3.));
        else if (deriv_y==1) {
            if (y < 1./3.)      return   0.5*std::exp(-0.5*fabs(y-1./3.));
            else                return  -0.5*std::exp(-0.5*fabs(y-1./3.));
        }
        else                    return   0.25*std::exp(-0.5*fabs(y-1./3.));
    }
    else if (nr==3) {
        if (deriv_y==0)         return  std::exp(-0.5*(y-0.1)*(y-0.1));
        else if (deriv_y==1)    return -std::exp(-0.5*(y-0.1)*(y-0.1)) * 2*0.5*(y-0.1);
        else                    return  std::exp(-0.5*(y-0.1)*(y-0.1))
                                        * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
    }
    else {
        assert(0); return 0;
    }
    return 0;
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::exact_z(T z, int deriv_z)
{
    if (nr==1) {
        if (deriv_z==0)         return  std::exp(-(z-0.1)*(z-0.1));
        else if (deriv_z==1)    return -std::exp(-(z-0.1)*(z-0.1)) * 2*(z-0.1);
        else                    return  std::exp(-(z-0.1)*(z-0.1)) * (4*(z-0.1)*(z-0.1)-2);
    }
    else if (nr==2) {
        if (deriv_z==0)         return  std::exp(-0.1*fabs(z-1./3.));
        else if (deriv_z==1) {
            if (z < 1./3.)      return   0.1*std::exp(-0.1*fabs(z-1./3.));
            else                return  -0.1*std::exp(-0.1*fabs(z-1./3.));
        }
        else                    return   0.01*std::exp(-0.1*fabs(z-1./3.));
    }
    else if (nr==3) {
        if (deriv_z==0)         return  std::exp(-2.*fabs(z-1./3.));
        else if (deriv_z==1) {
            if (z < 1./3.)      return   2.*std::exp(-2.*fabs(z-1./3.));
            else                return  -2.*std::exp(-2.*fabs(z-1./3.));
        }
        else                    return   4.*std::exp(-2.*fabs(z-1./3.));
    }
    else {
        assert(0); return 0;
    }
    return 0;
}

template <typename T>
T
TensorRefSols_PDE_Realline3D<T>::H1norm()
{
    T ret = 0.;
    if (nr==1) {
        T i11   = 3.963327297606011;
        T i22   = 1.772453850905516;
        T i33   = 1.2533141373155;
        T ddi11 = 0.3963327297606013;
        T ddi22 = 0.886226925452758;
        T ddi33 = 1.2533141373155;

        ret = sqrt(ddi11*i22*i33 + i11*ddi22*i33 + i11*i22*ddi33 + i11*i22*i33);
    }
    else if (nr==2) {
        T i11   = 0.5;
        T i22   = 2.;
        T i33   = 10.;
        T ddi11 = 2.;
        T ddi22 = 0.5;
        T ddi33 = 0.1;

        ret = sqrt(ddi11*i22*i33 + i11*ddi22*i33 + i11*i22*ddi33 + i11*i22*i33);
    }
    else if (nr==3) {
        T i11   = 3.963327297606011;
        T i22   = 1.772453850905516;
        T i33   = 0.5;
        T ddi11 = 0.3963327297606013;
        T ddi22 = 0.886226925452758;
        T ddi33 = 2.;

        ret = sqrt(ddi11*i22*i33 + i11*ddi22*i33 + i11*i22*ddi33 + i11*i22*i33);
    }
    return ret;
}

}   //namespace lawa
