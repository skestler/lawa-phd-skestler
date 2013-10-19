namespace lawa {

template <typename T>
int
RefSols_PDE_Realline2D<T>::nr;

template <typename T>
T
RefSols_PDE_Realline2D<T>::c;

template <typename T>
flens::DenseVector<Array<T> >
RefSols_PDE_Realline2D<T>::sing_pts_x;

template <typename T>
flens::DenseVector<Array<T> >
RefSols_PDE_Realline2D<T>::sing_pts_y;

template <typename T>
void
RefSols_PDE_Realline2D<T>::setExample(int _nr, T _c)
{
    c=_c;
    assert(c>=0);
    nr=_nr;

    if (nr==1) {
        sing_pts_x.engine().resize(3);
        sing_pts_y.engine().resize(3);
        sing_pts_x = -1., 0., 1.;
        sing_pts_y = -1., 0., 1.;
    }
    if (nr==2) {
/*
        sing_pts_x.engine().resize(1);
        sing_pts_y.engine().resize(1);
        sing_pts_x = 0.1;
        sing_pts_y = 0.1;
*/
        sing_pts_x.engine().resize(9);
        sing_pts_y.engine().resize(9);
        sing_pts_x = -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5;
        sing_pts_y = -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5;

/*
        sing_pts_x.engine().resize(17);
        sing_pts_y.engine().resize(17);
        sing_pts_x = -0.3,-0.2, -0.1, 0., 0.099, 0.0999,  0.09999, 0.09999999, 0.1, 0.10000001, 0.10001, 0.1001, 0.101, 0.2, 0.3, 0.4, 0.5;
        sing_pts_y = -0.3,-0.2, -0.1, 0., 0.099, 0.0999,  0.09999, 0.09999999, 0.1, 0.10000001, 0.10001, 0.1001, 0.101, 0.2, 0.3, 0.4, 0.5;
*/
    }
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::exact(T x, T y)
{
    return exact(x,y,0,0);
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::minus_exact(T x, T y)
{
    return -exact(x,y,0,0);
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::exact_dx(T x, T y)
{
    return exact(x,y,1,0);
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::exact_dy(T x, T y)
{
    return exact(x,y,0,1);
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::rhs(T x, T y)
{
    return -exact(x,y,2,0)-exact(x,y,0,2)+c*exact(x,y,0,0);
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::exact(T x, T y, int deriv_x, int deriv_y)
{
    if (nr==1) {
        if ((deriv_x==0) && (deriv_y==0)) {
            return exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
        }
        else if ((deriv_x==0) && (deriv_y==2)) {
            return (pow(-(x-0.1)-2*(y-0.1),2) - 2)*
                   exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
        }
        else if ((deriv_x==2) && (deriv_y==0)) {
            return (pow(-4*(x-0.1)-(y-0.1),2) - 4)*
                   exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
        }
        else {
            return 0.;
        }
    }
    else if (nr==2) {
        T XmA_p2_p_YmB_p2 = (x-0.1)*(x-0.1)+(y-0.1)*(y-0.1);
        if ((deriv_x==0) && (deriv_y==0)) {
            return exp(-sqrt( XmA_p2_p_YmB_p2 ));
        }
        else if ((deriv_x==0) && (deriv_y==1)) {
            return (-(y-0.1)*sqrt(XmA_p2_p_YmB_p2)/XmA_p2_p_YmB_p2)*exp(-sqrt( XmA_p2_p_YmB_p2 ));
        }
        else if ((deriv_x==1) && (deriv_y==0)) {
            return (-(x-0.1)*sqrt(XmA_p2_p_YmB_p2)/XmA_p2_p_YmB_p2)*exp(-sqrt( XmA_p2_p_YmB_p2 ));
        }
        else if ((deriv_x==0) && (deriv_y==2)) {
            return ( (y-0.1)*(y-0.1)/XmA_p2_p_YmB_p2 -
                     1./sqrt(XmA_p2_p_YmB_p2) + (y-0.1)*(y-0.1)/pow(XmA_p2_p_YmB_p2,1.5) )*
                   exp(-sqrt( XmA_p2_p_YmB_p2 ));
        }
        else if ((deriv_x==2) && (deriv_y==0)) {
            return ( (x-0.1)*(x-0.1)/XmA_p2_p_YmB_p2 -
                     1./sqrt(XmA_p2_p_YmB_p2) + (x-0.1)*(x-0.1)/pow(XmA_p2_p_YmB_p2,1.5) )*
                   exp(-sqrt( XmA_p2_p_YmB_p2 ));
        }
        else {
            return 0.;
        }
    }
    else {
        assert(0);
        return 0;
    }
}

template <typename T>
T
RefSols_PDE_Realline2D<T>::H1norm()
{
    T ret = 0.;
    if (nr==1)             {
        ret = 1.187410411723726 + 2.374820823447452 + 1.187410411723726;
        ret = sqrt(ret);
    }
    else if (nr==2) {
        //ret = 1.570795505591071 + 0.7853981678060725 + 0.7853981678060725;
        //ret = 0.5*M_PI + 0.7853982269 + 0.7853982269;
        ret = M_PI;
        ret = sqrt(ret);
    }
    return ret;
}

}   // namespace lawa
