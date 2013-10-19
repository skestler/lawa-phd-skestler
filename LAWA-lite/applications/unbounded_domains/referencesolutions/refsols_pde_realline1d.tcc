namespace lawa {

template <typename T>
int
RefSols_PDE_Realline1D<T>::nr;

template <typename T>
T
RefSols_PDE_Realline1D<T>::reaction;

template <typename T>
T
RefSols_PDE_Realline1D<T>::convection;

template <typename T>
T
RefSols_PDE_Realline1D<T>::diffusion;

template <typename T>
flens::DenseVector<Array<T> >
RefSols_PDE_Realline1D<T>::sing_pts;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
RefSols_PDE_Realline1D<T>::deltas;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
RefSols_PDE_Realline1D<T>::H1_deltas;


template <typename T>
void
RefSols_PDE_Realline1D<T>::setExample(int _nr, T _reaction, T _convection, T _diffusion)
{
    nr=_nr;
    reaction      = _reaction;
    convection    = _convection;
    diffusion     = _diffusion;

    assert(nr>=1);
    assert(nr<=6);
    if (nr==1) {

    }
    else if (nr==2) {
        sing_pts.engine().resize(1);
        sing_pts(1) = 0.01;
        deltas.engine().resize(1,2);
        deltas(1,1) = 0.01; deltas(1,2) = diffusion*2.;
        H1_deltas.engine().resize(1,2);
        H1_deltas(1,1) = 0.01; H1_deltas(1,2) = 2.;
    }
    else if (nr==3) {
        sing_pts.engine().resize(2);
        sing_pts(1) = -M_PI/32.;
        sing_pts(2) =  M_PI/32.;
        deltas.engine().resize(2,2);
        //deltas(1,1) = -M_PI/32.; deltas(1,2) = -diffusion*4.;
        //deltas(2,1) =  M_PI/32.; deltas(2,2) = diffusion*4.;
        deltas(1,1) = -M_PI/32.; deltas(1,2) = -diffusion*0.4;
        deltas(2,1) =  M_PI/32.; deltas(2,2) = diffusion*0.1;
        H1_deltas.engine().resize(2,2);
        H1_deltas(1,1) = -M_PI/32.; H1_deltas(1,2) = -0.4;
        H1_deltas(2,1) =  M_PI/32.; H1_deltas(2,2) =  0.1;
    }
    else if (nr==4) {    //second derivative not continuous on R!!
        sing_pts.engine().resize(2);
        sing_pts(1) = -0.4;
        sing_pts(2) =  0.9;
    }
    else if (nr==5) {    //second derivative not continuous on R!!
        sing_pts.engine().resize(2);
        sing_pts(1) = -0.4;
        sing_pts(2) =  0.9;
        deltas.engine().resize(2,2);
        deltas(1,1) = -0.4; deltas(1,2) = diffusion*20.8;
        deltas(2,1) =  0.9; deltas(2,2) = diffusion*105.3;
        H1_deltas.engine().resize(2,2);
        H1_deltas(1,1) = -0.4; H1_deltas(1,2) = 20.8;
        H1_deltas(2,1) =  0.9; H1_deltas(2,2) = 105.3;
    }
    else if (nr==6) {
        sing_pts.engine().resize(1);
        sing_pts(1) = 0.0001;
        deltas.engine().resize(1,2);
        //deltas(1,1) = 0.0001; deltas(1,2) = 3.;
        deltas(1,1) = 0.0001; deltas(1,2) = diffusion*9.;
        H1_deltas.engine().resize(1,2);
        H1_deltas(1,1) = 0.0001; H1_deltas(1,2) = 9.;
    }
    else if (nr==7) {
        sing_pts.engine().resize(1);
        sing_pts(1) = 0.0001;
        deltas.engine().resize(1,2);
        deltas(1,1) = 0.0001; deltas(1,2) = diffusion*3.;
        H1_deltas.engine().resize(1,2);
        H1_deltas(1,1) = 0.0001; H1_deltas(1,2) = 3.;
    }
    else if (nr==8) {
        sing_pts.engine().resize(2);
        sing_pts(1) = -16./3.;
        sing_pts(2) =  16./3.;
        deltas.engine().resize(2,2);
        deltas(1,1) = -16./3.; deltas(1,2) = diffusion*22773./6391.;
        deltas(2,1) =  16./3.; deltas(2,2) = diffusion*16382./6391.;
        H1_deltas.engine().resize(2,2);
        H1_deltas(1,1) = -16./3.; H1_deltas(1,2) = 22773./6391.;
        H1_deltas(2,1) =  16./3.; H1_deltas(2,2) = 16382./6391.;
    }
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::exact(T x, int deriv)
{

    if      (nr==1) {
        if (deriv == 0)         return 10.*std::exp(-0.1*(x-0.1)*(x-0.1));
        else if (deriv == 1)    return -2*(x-0.1)*std::exp(-0.1*(x-0.1)*(x-0.1));
        else if (deriv == 2)    return (4*0.1*(x-0.1)*(x-0.1)-2)*std::exp(-0.1*(x-0.1)*(x-0.1));
        else                    assert(0);
    }
    else if (nr==2) {
        if (deriv == 0)             return 10.*std::exp(-0.1*fabs(x-0.01));
        else if (deriv == 1) {
            if (x > 0.01)         return -std::exp(-0.1*(x-0.01));
            else                  return  std::exp( 0.1*(x-0.01));
        }
        else if (deriv == 2) {    //only piecewise, delta is required!!!
            if (x > 0.01)         return 0.1*std::exp(-0.1*(x-0.01));
            else                  return 0.1*std::exp( 0.1*(x-0.01));
        }
        else                    assert(0);
    }
    else if (nr == 3){
        if (deriv == 0) {
            if (x<-M_PI/32.)          return -1./std::pow(0.1*(x+M_PI/32.)-1,(T)4);
            else if (x>M_PI/32.)      //return  1./std::pow(x-M_PI/32.+1,(T)4);
                                    return  std::exp(-0.1*(x-M_PI/32));
            else                     return std::sin(16*x);
         }
         else if (deriv == 1) {
             if (x<-M_PI/32.)         return  0.4/std::pow(0.1*(x+M_PI/32.)-1,(T)5);
             else if (x>M_PI/32.)     //return -4./std::pow(x-M_PI/32.+1,(T)5);
                                     return -0.1*std::exp(-0.1*(x-M_PI/32));
             else                    return 16*std::cos(16*x);
         }
         else if (deriv == 2) {
             if (x<-M_PI/32.)         return -0.2/std::pow(0.1*(x+M_PI/32.)-1,(T)6);
             else if (x>M_PI/32.)     //return  20./std::pow(x-M_PI/32.+1,(T)6);
                                     return 0.01*std::exp(-0.1*(x-M_PI/32));
             else                    return -256*std::sin(16*x);

         }
         else                         assert(0);
    }
    else if (nr==4) {
        if (deriv==0) {
            if ((x>-0.4) && (x<0.9))  return 100*(x+0.4)*(x+0.4)*(x-0.9)*(x-0.9);
            else                       return 0;
        }
        else if (deriv==1) {
            if ((x>-0.4) && (x<0.9))  return   200*(x+0.4)*(x-0.9)*(x-0.9)
                                             + 200*(x+0.4)*(x+0.4)*(x-0.9);
            else                       return 0;
        }
        else if (deriv==2) {
            if ((x>-0.4) && (x<0.9))  return   200*(x-0.9)*(x-0.9) + 800.*(x+0.4)*(x-0.9)
                                             + 200.*(x+0.4)*(x+0.4);
            else                       return 0;
        }
        else                         assert(0);
    }
    else if (nr==5) {
        if (deriv==0) {
            if ((x>-0.4) && (x<0.9))   return 100*(x+0.4)*x*x*(x-0.9);
            else return 0;
        }
        else if (deriv==1) {
            if ((x>-0.4) && (x<0.9))   return   100*(x-0.9)*x*x + 200*(x-0.9)*x*(x+0.4)
                                              + 100*x*x*(x+0.4);
            else return 0;
        }
        else if (deriv==2) {
            if ((x>-0.4) && (x<0.9))   return   400*x*(x-0.9) + 200*x*x + 200*(x-0.9)*(x+0.4)
                                              + 400*(x+0.4)*x;
            else                       return 0;
        }
        else                         assert(0);
    }
    else if (nr==6) {
        if (deriv==0) {
            if (x<=0.0001)      return 1./(std::pow(x-0.0001-1,(T)4));
            else                return 1./(std::pow(x-0.0001+1,(T)5));
        }
        else if (deriv==1) {
            if (x<=0.0001)      return -4./(std::pow(x-0.0001-1,(T)5));
            else                return -5./(std::pow(x-0.0001+1,(T)6));
        }
        else if (deriv==2) {
            if (x<=0.0001)      return 20./(std::pow(x-0.0001-1,(T)6));
            else                return 30./(std::pow(x-0.0001+1,(T)7));
        }
    }
    else if (nr==7) {
        T tmp = x-0.0001-1;
        if (deriv==0) {
            if (x<=0.0001)      return exp(2*(x-0.0001));
            else                return exp(-(x-0.0001));
        }
        else if (deriv==1) {
            if (x<=0.0001)      return 2*exp(2*(x-0.0001));
            else                return -exp(-(x-0.0001));
        }
        else if (deriv==2) {
            if (x<=0.0001)      return 4*exp(2*(x-0.0001));
            else                return exp(-(x-0.0001));
        }
    }
    else if (nr == 8){
        if (deriv == 0) {
            T coeff1 = 675./12782., coeff2 = -3209./6391.;
            if (x<-16./3.)         return -1./std::pow(x-1+16./3.,(T)3.);
            else if (x>16/3.)      return  1./std::pow(x+1-16./3.,(T)2.);
            else                 return coeff1*x*x+coeff2;
        }
        else if (deriv == 1) {
            T coeff1 = 675./12782.;
            if (x<-16./3.)             return  3./std::pow(x-1+16./3.,(T)4.);
            else if (x>16./3.)      return -2./std::pow(x+1-16./3.,(T)3.);
            else                     return 2.*coeff1*x;
        }
        else if (deriv == 2) {
            T coeff1 = 675./12782.;
            if (x<-16./3.)             return  -12./std::pow(x-1+16./3.,(T)5.);
            else if (x>16./3.)      return 6./std::pow(x+1-16./3.,(T)4.);
            else                     return 2.*coeff1;
        }
        else                         assert(0);
    }
    else {
        assert(0);
        return 0;
    }
    assert(0);
    return 0;
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::u(T x)
{
    return RefSols_PDE_Realline1D<T>::exact(x, 0);
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::d_u(T x)
{
    return RefSols_PDE_Realline1D<T>::exact(x, 1);
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::rhs(T x)
{
    return -diffusion*exact(x,2) + convection*exact(x,1) + reaction*exact(x,0);
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::H1_rhs(T x)
{
    return -exact(x,2) + exact(x,0);
}

template <typename T>
T
RefSols_PDE_Realline1D<T>::H1norm()
{
    if (nr==1)          return   std::sqrt(100*std::sqrt(0.5*M_PI/(0.1))
                               + std::sqrt(0.5*M_PI)/std::pow(0.1,(T)1.5) );
    else if (nr==2)     return std::sqrt(1000. + 10.);
    //else if (nr==3)     return std::sqrt( (16+7.*M_PI)/56. + (2./9.)*(16.+9.*M_PI)   );
    //else if (nr==3)     return std::sqrt( (64.+7.*M_PI)/224. + (8./9.)*(4.+9.*M_PI)   );
    //else if (nr==3)     return std::sqrt( (1152.+7.*M_PI)/224. + (1./180.)*(329.+1440.*M_PI)   );
    else if (nr==3)     return std::sqrt( (1440.+7.*M_PI)/224. + (1./180.)*(41.+1440.*M_PI)   );
    else if (nr==4)     return std::sqrt(168.3253868730168 + 1195.209847619049);
    else if (nr==5)     return std::sqrt(43.3676117539685 + 980.929114285711);
    else if (nr==6)     return std::sqrt(0.253968253968254 + 4.05050505050505);
    else if (nr==7)     return std::sqrt(0.75 + 1.5);
    else if (nr==8)     return std::sqrt(1637492008./612673215.  + 656353759./204224405.);
    else {
        assert(0);
        return 0;
    }
    assert(0);
    return 0;

}

template <typename T>
int
RefSols_PDE_Realline1D<T>::getMinimalLevel(int d, int d_)
{
    if (diffusion!=1 || convection!=0 || reaction!=1) {
        std::cerr << "Right-hand side parameters only implemented for -u'' + u" << std::endl;
        assert(0);
    }
    if (d==2 && d_==2) {
        if (nr==1)       return -2;
        else if (nr==2)  return -4;
        else if (nr==3)  return -1;
        else if (nr==4)  return  0;
        else if (nr==5)  return  0;
        else if (nr==6)  return -1;//better convergence behaviour, works also for higher levels
        else {  assert(0); return 0; }
    }
    else if (d==3 && d_==3) {
        if (nr==1)       return -2;
        else if (nr==2)  return -4;
        else if (nr==3)  return -1;
        else if (nr==4)  return  0;
        else if (nr==5)  return -1;
        else if (nr==6)  return -1;//better convergence behaviour, works also for higher levels
        else {  assert(0); return 0; }
    }
    else if (d==3 && d_==5) {
        if (nr==1)       return -2;
        else if (nr==2)  return -4;
        else if (nr==3)  return -1;
        else if (nr==4)  return  0;
        else if (nr==5)  return -1;
        else if (nr==6)  return -1;//better convergence behaviour, works also for higher levels
        else {  assert(0); return 0; }
    }
    else {
        assert(0);
        return 0;
    }
    assert(0);
    return 0;
}

template <typename T>
void
RefSols_PDE_Realline1D<T>::getRHS_W_XBSplineParameters(int d, int d_,
                                                       T &_left_bound, T &_right_bound,
                                                       int &_J_plus_smooth, int &_J_plus_singular,
                                                       bool &_singular_integral, T /*eps*/)
{
    if (diffusion!=1 || convection!=0 || reaction!=1) {
        std::cerr << "Right-hand side parameters only implemented for -u'' + u" << std::endl;
        assert(0);
    }
    if (nr == 1) {
        _singular_integral=false;
        if (d==2 && d_==2) {
            _left_bound = -20.; _right_bound = 20.;
            _J_plus_smooth = 7;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _left_bound = -80.; _right_bound = 80.;
            _J_plus_smooth = 7;
            _J_plus_singular = 50;
        }
        else if (d==4) {
            _left_bound = -40.; _right_bound = 40.;
            _J_plus_smooth = 7;
            _J_plus_singular = 50;
        }
    }
    else if (nr == 2) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -114.; _right_bound = 114.;
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _left_bound = -214.; _right_bound = 214.;
            _J_plus_smooth = 5;
            _J_plus_singular = 40;
        }
        else if (d==4) {
            _left_bound = -220.; _right_bound = 220.;
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
    }
    else if (nr == 3) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -100.; _right_bound = 100.;
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _left_bound = -100.; _right_bound = 100.;
            _J_plus_smooth = 4;
            _J_plus_singular = 60;
        }
        else if (d==4) {
            _left_bound = -120.; _right_bound = 120.;
            _J_plus_smooth = 5;
            _J_plus_singular = 60;
        }
    }
    else if (nr == 4) {
        _singular_integral=true;
        _left_bound = -0.4; _right_bound = 0.9;
        if (d==2 && d_==2) {
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _J_plus_smooth = 5;
            _J_plus_singular = 50;
        }
        else if (d==4) {
            _J_plus_smooth = 5;
            _J_plus_singular = 50;
        }
    }
    else if (nr == 5) {
        _singular_integral=true;
        _left_bound = -0.4; _right_bound = 0.9;
        if (d==2 && d_==2) {
            _J_plus_smooth = 8;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
        else if (d==4) {
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
    }
    else if (nr == 6) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -114.; _right_bound = 114.;
            _J_plus_smooth = 6;
            _J_plus_singular = 40;
        }
        else if (d==3) {
            _left_bound = -200.; _right_bound = 200.;
            _J_plus_smooth = 3;
            _J_plus_singular = 40;
        }
        else if (d==4) {
            _left_bound = -300.; _right_bound = 300.;
            _J_plus_smooth = 5;
            _J_plus_singular = 40;
        }
    }
    else {
        std::cerr << "Parameters not set for example " << nr << " and d=" << d
                  << ", d_=" << d_ << std::endl;
        assert(0);
    }

}


template <typename T>
void
RefSols_PDE_Realline1D<T>::getRHS_WO_XBSplineParameters(int d, int d_,
                                                        T &_left_bound, T &_right_bound,
                                                        int &_J_plus_smooth, int &_J_minus_smooth,
                                                        int &_J_plus_singular,
                                                        int &_J_minus_singular,
                                                        bool &_singular_integral, T /*eps*/)
{
    if (diffusion!=1 || convection!=0 || reaction!=1) {
        std::cerr << "Right-hand side parameters only implemented for -u'' + u" << std::endl;
        assert(0);
    }
    if (nr == 1) {
        _singular_integral=false;
        if (d==2 && d_==2) {
            _left_bound = -40.; _right_bound = 40.;
            _J_plus_smooth = 6;     _J_minus_smooth = -40;
            _J_plus_singular = 40;     _J_minus_singular = -40;
        }
        else if (d==3) {
            _left_bound = -60.; _right_bound = 60.;
            _J_plus_smooth = 4;     _J_minus_smooth = -50;
            _J_plus_singular = 50;  _J_minus_singular = -50;
        }
    }
    else if (nr == 2) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -400.; _right_bound = 400.;
            _J_plus_smooth = 7;     _J_minus_smooth = -50;
            _J_plus_singular = 40;     _J_minus_singular = -50;
        }
        else if (d==3) {
            _left_bound = -200.; _right_bound = 200.;
            _J_plus_smooth = 4;     _J_minus_smooth = -80;
            _J_plus_singular = 80;     _J_minus_singular = -80;
        }
    }
    else if (nr == 3) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -300.; _right_bound = 300.;
            _J_plus_smooth = 6;     _J_minus_smooth = -40;
            _J_plus_singular = 40;     _J_minus_singular = -40;
        }
        else if (d==3) {
            _left_bound = -300.; _right_bound = 300.;
            _J_plus_smooth = 5;     _J_minus_smooth = -80;
            _J_plus_singular = 40;     _J_minus_singular = -80;
        }
    }
    else if (nr == 4) {
        _singular_integral=true;
        _left_bound = -0.4; _right_bound = 0.9;
        if (d==2 && d_==2) {
            _J_plus_smooth = 6;     _J_minus_smooth = -40;
            _J_plus_singular = 40;     _J_minus_singular = -40;
        }
        else if (d==3) {
            _J_plus_smooth = 5;     _J_minus_smooth = -50;
            _J_plus_singular = 50;     _J_minus_singular = -50;
        }
    }
    else if (nr == 5) {
        _singular_integral=true;
        _left_bound = -0.4; _right_bound = 0.9;
        if (d==2 && d_==2) {
            _J_plus_smooth = 10;     _J_minus_smooth = -40;
            _J_plus_singular = 40;     _J_minus_singular = -40;
        }
        else if (d==3) {
            _J_plus_smooth = 11;     _J_minus_smooth = -50;
            _J_plus_singular = 55;     _J_minus_singular = -50;
        }
    }
    else if (nr == 6) {
        _singular_integral=true;
        if (d==2 && d_==2) {
            _left_bound = -114.; _right_bound = 114.;
            _J_plus_smooth = 6;     _J_minus_smooth = -40;
            _J_plus_singular = 45;     _J_minus_singular = -40;
        }
        else if (d==3) {
            _left_bound = -200.; _right_bound = 200.;
            _J_plus_smooth = 5;     _J_minus_smooth = -50;
            _J_plus_singular = 46;     _J_minus_singular = -50;
        }
    }
    else {
        std::cerr << "Parameters not set for example " << nr << " and d=" << d
                  << ", d_=" << d_ << std::endl;
        assert(0);
    }
}

} //namespace lawa
