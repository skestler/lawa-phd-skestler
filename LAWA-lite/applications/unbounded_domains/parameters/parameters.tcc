namespace lawa {

template <typename T>
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1NormPreconditioner1D<T, Basis<T,Primal,R,CDF> > >::Parameters
           (const Basis<T,Primal,R,CDF> &_basis,
            const HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > &_Bil, bool _w_XBSpline, int _j0)
           : basis(_basis), Bil(_Bil), w_XBSpline(_w_XBSpline), j0(_j0),
             cA(0.), CA(0.), kappa(0.),
             alpha(0.), omega(0.), gamma(0.), theta(0.)
{
    if (!_w_XBSpline) {    //only wavelet discretization
         std::cout << "Parameters: j0 = - Infinity" << std::endl;
        if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
            cA = 0.19; CA = 2.;
        }
        else if (basis.d==3 && basis.d_==3 && Bil.c == 1.) {
            cA = 0.038; CA = 2.7;
        }
        else if (basis.d==3 && basis.d_==5 && Bil.c == 1.) {
            cA = 0.16; CA = 2.8;
        }
        else assert(0);
    }
    else {                //wavelet with b-splines on coarse level discretization
         std::cout << "Parameters: j0 = " << j0 << std::endl;
        if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
            if (j0==0)  {    cA = 0.375; CA = 2.1;    }
            else if (j0==-1) {    cA = 0.58;  CA = 1.86;    }
            else if (j0==-2) {    cA = 0.58;  CA = 1.86;    }
            else if (j0==-3) {    cA = 0.55;  CA = 1.86;    }
            else if (j0==-4) {    cA = 0.46;  CA = 1.86;    }
            else if (j0==-5) {    cA = 0.4;   CA = 1.86;    }
            else if (j0==-6) {    cA = 0.36;  CA = 1.89;    }
            else if (j0==-7) {    cA = 0.33;  CA = 1.92;    }
            else if (j0==-8) {    cA = 0.3;   CA = 1.95;    }
            else if (j0==-8) {    cA = 0.29;  CA = 1.95;    }
            else if (j0==-10) {    cA = 0.27;  CA = 1.96;    }
            else assert(0);
        }
        else if (basis.d==3 && basis.d_==3 && Bil.c == 1.) {
            if (j0==0)  {    cA = 0.43;  CA = 1.94;    }
            else if (j0==-1) {    cA = 0.39;  CA = 2.03;    }
            else if (j0==-2) {    cA = 0.29;  CA = 2.24;    }
            else if (j0==-3) {    cA = 0.21;  CA = 2.44;    }
            else if (j0==-4) {    cA = 0.15;  CA = 2.55;    }
            else if (j0==-5) {    cA = 0.12;  CA = 2.59;    }
            else if (j0==-6) {    cA = 0.11;   CA = 2.61;    }
            else if (j0==-7) {    cA = 0.09;  CA = 2.62;    }
            else if (j0==-8) {    cA = 0.078; CA = 2.63;    }
            else if (j0==-9) {    cA = 0.07;  CA = 2.64;    }
            else if (j0==-10) {    cA = 0.063; CA = 2.64;    }
            else assert(0);
        }
        else if (basis.d==3 && basis.d_==5 && Bil.c == 1.) {
            if (j0==0)  {    cA = 0.457;  CA = 1.96;    }
            else if (j0==-1) {    cA = 0.416;  CA = 2.07;    }
            else if (j0==-2) {    cA = 0.326;  CA = 2.32;    }
            else if (j0==-3) {    cA = 0.24;   CA = 2.53;    }
            else if (j0==-4) {    cA = 0.197;  CA = 2.66;    }
            else if (j0==-5) {    cA = 0.18;   CA = 2.7;    }
            else if (j0==-6) {    cA = 0.17;   CA = 2.71;    }
            else assert(0);
        }
        else assert(0);
    }
    kappa = CA/cA;
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

    assert(cA>0); assert(CA>0); assert(omega>0.);
    assert(alpha>0.); assert(gamma>0.); assert(theta>0.);
    assert((alpha+omega)/(1-omega)<1./std::sqrt(kappa));
}

template <typename T>
void
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1NormPreconditioner1D<T, Basis<T,Primal,R,CDF> > >::getGHSADWAVParameters
                                                                (T &_alpha, T &_omega,
                                                                 T &_gamma, T &_theta) const
{
    _alpha = alpha; _omega = omega; _gamma = gamma; _theta = theta;
}


template <typename T>
Parameters<T, Basis<T,Orthogonal,R,Multi>, HelmholtzOperator1D<T,Basis<T,Orthogonal,R,Multi> >,
           H1NormPreconditioner1D<T, Basis<T,Orthogonal,R,Multi> > >::Parameters
           (const Basis<T,Orthogonal,R,Multi> &_basis,
            const HelmholtzOperator1D<T,Basis<T,Orthogonal,R,Multi> > &_Bil,int _j0)
           : basis(_basis), Bil(_Bil), w_XBSpline(true), j0(_j0),
             cA(0.), CA(0.), kappa(0.),
             alpha(0.), omega(0.), gamma(0.), theta(0.)
{

    if (basis.d==2 && Bil.c == 1.) {
        if (j0==-1)      {    cA = 0.14;  CA = 2.9;    }
        else if (j0==-2) {    cA = 0.19;  CA = 2.9;    }
        else if (j0==-3) {    cA = 0.19;  CA = 2.9;    }
        else if (j0==-4) {    cA = 0.19;  CA = 2.9;    }
        else assert(0);
    }
    else assert(0);

    kappa = CA/cA;
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

    assert(cA>0); assert(CA>0); assert(omega>0.);
    assert(alpha>0.); assert(gamma>0.); assert(theta>0.);
    assert((alpha+omega)/(1-omega)<1./std::sqrt(kappa));
}

template <typename T>
void
Parameters<T, Basis<T,Orthogonal,R,Multi>, HelmholtzOperator1D<T,Basis<T,Orthogonal,R,Multi> >,
           H1NormPreconditioner1D<T, Basis<T,Orthogonal,R,Multi> > >::getGHSADWAVParameters
                                                                      (T &_alpha, T &_omega,
                                                                       T &_gamma, T &_theta) const
{
    _alpha = alpha; _omega = omega; _gamma = gamma; _theta = theta;
}


template <typename T>
Parameters<T, Basis<T,Primal,R,CDF>, WeightedHelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           WeightedSobolevMidPointPreconditioner1D<T, Basis<T,Primal,R,CDF> > >::Parameters
           (const Basis<T,Primal,R,CDF> &_basis,
            const WeightedHelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > &_Bil,
            bool _w_XBSpline, int _j0)
           : basis(_basis), Bil(_Bil), w_XBSpline(_w_XBSpline), j0(_j0),
             cA(0.), CA(0.), kappa(0.),
             alpha(0.), omega(0.), gamma(0.), theta(0.)
{
    if (!_w_XBSpline) {    //only wavelet discretization
        assert(0);
    }
    else {                //wavelet with b-splines on coarse level discretization
         std::cout << "Parameters: j0 = " << j0 << std::endl;
        if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
            if      (j0==0)  {    cA = 0.064; CA = 8.9;    }
            else if (j0==1)  {    cA = 0.18; CA = 3.18;    }
            else if (j0==2)  {    cA = 0.24; CA = 2.28;    }
            else assert(0);
        }
        else if (basis.d==3 && basis.d_==3 && Bil.c == 1.) {
            assert(0);
        }
        else if (basis.d==3 && basis.d_==5 && Bil.c == 1.) {
            assert(0);
        }
        else assert(0);
    }
    kappa = CA/cA;
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

    assert(cA>0); assert(CA>0); assert(omega>0.);
    assert(alpha>0.); assert(gamma>0.); assert(theta>0.);
    assert((alpha+omega)/(1-omega)<1./std::sqrt(kappa));
}

template <typename T>
void
Parameters<T, Basis<T,Primal,R,CDF>, WeightedHelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           WeightedSobolevMidPointPreconditioner1D<T, Basis<T,Primal,R,CDF> > >::
           getGHSADWAVParameters(T &_alpha, T &_omega, T &_gamma, T &_theta) const
{
    _alpha = alpha; _omega = omega; _gamma = gamma; _theta = theta;
}

}
