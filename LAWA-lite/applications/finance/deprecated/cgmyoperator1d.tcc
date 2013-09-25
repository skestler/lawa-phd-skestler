namespace lawa {

template <typename T, typename Basis>
CGMYOperator1D<T,Basis>::CGMYOperator1D(const Basis &_basis, T _diffusion, T _convection, T _reaction,
                                        T C, T G, T M, T Y, int _order, int _n, T _sigma)
    : basis(_basis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      phi(basis.mra), d_phi(basis.mra,1), psi(basis), d_psi(basis,1),
      integral_sfsf(phi, phi), d_integral_sfsf(phi, d_phi), dd_integral_sfsf(d_phi, d_phi),
      integral_sfw(phi, psi),  d_integral_sfw(phi, d_psi),  dd_integral_sfw(d_phi, d_psi),
      integral_wsf(psi, phi),  d_integral_wsf(psi, d_phi),  dd_integral_wsf(d_psi, d_phi),
      integral_ww(psi,psi),    d_integral_ww(psi,d_psi),    dd_integral_ww(d_psi,d_psi),
      cgmy(C,G,M,Y), cgmy_adjoint(C,M,G,Y),
      n(_n), sigma(_sigma), omega(0.), mu(0.5), order(_order)
{
    assert(diffusion>=0.);
    assert(reaction >=0.);
    assert(C>=0);
    T C_f = 1.;
    omega = C_f*(1-sigma)/(4.*sigma);
    assert(omega<1);

    int max_order = std::max(order,int(std::ceil(n*mu)));

    std::cout << "CGMYOperator1d: n = " << n << ", sigma = " << sigma << ", omega = " << omega
              << ", mu = " << mu  << ", max_order = " << max_order << std::endl;

    T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(max_order, max_order);
    _weights.engine().resize(max_order, max_order);
    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=max_order; ++k) {
        int     m = (k+1)/2;
        T xm = 0.5 * (x2+x1),
          xl = 0.5 * (x2-x1);
        for (int i=1; i<=m; ++i) {
            T z = cos(M_PI*(i-0.25)/(k+0.5)),
              z1, pp;
            do {
                T p1 = 1.0,
                  p2 = 2.0;
                for (int j=1; j<=k; ++j) {
                    T p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            _knots(k,i)     = xm - xl*z;
            _knots(k,k+1-i) = xm + xl*z;
            _weights(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
            _weights(k,k+1-i) = _weights(k,i);
        }
    }
}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    if (row_index.j >= col_index.j) {
        return CGMYOperator1D<T,Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                   col_index.xtype, col_index.j, col_index.k,
                                                   cgmy);
    }
    else {
        return CGMYOperator1D<T,Basis>::operator()(col_index.xtype, col_index.j, col_index.k,
                                                   row_index.xtype, row_index.j, row_index.k,
                                                   cgmy_adjoint);
    }
}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::getc() const
{
    return reaction;
}

template <typename T, typename Basis>
const Basis&
CGMYOperator1D<T,Basis>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::operator()(XType xtype_row, int j_row, int k_row,
                                    XType xtype_col, int j_col, int k_col,
                                    const CGMYUtils<T> &cgmy_comp) const
{
    //PDE part
    T val = 0;
    T d_val = 0;
    T dd_val = 0;

    if(xtype_row == XBSpline){
        if(xtype_col == XBSpline){
             val = integral_sfsf(j_row, k_row, j_col, k_col);
             d_val = d_integral_sfsf(j_row, k_row, j_col, k_col);
             dd_val = dd_integral_sfsf(j_row, k_row, j_col, k_col);
        }
        else{
             val = integral_sfw(j_row, k_row, j_col, k_col);
             d_val = d_integral_sfw(j_row, k_row, j_col, k_col);
             dd_val = dd_integral_sfw(j_row, k_row, j_col, k_col);
        }
    }
    else{
        if(xtype_col == XBSpline){
             val = integral_wsf(j_row, k_row, j_col, k_col);
             d_val = d_integral_wsf(j_row, k_row, j_col, k_col);
             dd_val = dd_integral_wsf(j_row, k_row, j_col, k_col);
        }
        else{
            val = integral_ww(j_row, k_row, j_col, k_col);
            d_val = d_integral_ww(j_row, k_row, j_col, k_col);
            dd_val = dd_integral_ww(j_row, k_row, j_col, k_col);
        }
    }

    //Integral part
    T int_val = 0.;
    if (cgmy.C == 0) return diffusion*dd_val + convection * d_val + reaction*val + int_val;

    flens::DenseVector<Array<T> > singularPoints;
    if (xtype_row == XBSpline)  singularPoints = phi.singularSupport(j_row,k_row);
    else                         singularPoints = psi.singularSupport(j_row,k_row);

    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        T a = singularPoints(i);
        T b = singularPoints(i+1);
        if (a==b) continue;

        int_val += _integrate_dpsi_vs_int_dpsi_k(a,  b, xtype_row, j_row, k_row,
                                                        xtype_col, j_col, k_col, cgmy_comp);
    }
    return diffusion*dd_val + convection * d_val + reaction*val + int_val;
}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::_integrate_dpsi_vs_int_dpsi_k(T a, T b, XType xtype_row, int j_row, int k_row,
                                                                 XType xtype_col, int j_col, int k_col,
                                                                 const CGMYUtils<T> &cgmy_comp) const
{
    T ret = 0.0;
    flens::DenseVector<Array<T> > singularPoints_psi_col;
    if (xtype_col == XBSpline) {    singularPoints_psi_col = phi.singularSupport(j_col,k_col);    }
    else {                            singularPoints_psi_col = psi.singularSupport(j_col,k_col);    }

    for (int i=1; i<=order; ++i) {
        T x_ast = 0.5*(b-a)*_knots(order,i)+0.5*(b+a);
        T tmp = 0.;

        for (int j=singularPoints_psi_col.firstIndex(); j<singularPoints_psi_col.lastIndex(); ++j) {
            T a_y = singularPoints_psi_col(j)  -x_ast;
            T b_y = singularPoints_psi_col(j+1)-x_ast;
            //assert(a_y!=0.); assert(b_y!=0.);
            if ((a_y < 0.) && (b_y >0.)) {
                if (cgmy.Y<1) {
                    tmp += _nonsingular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, 0., xtype_col, j_col, k_col, cgmy_comp);
                    tmp += _nonsingular_quadrature_dpsi_vs_CGMYkernel(x_ast, 0., b_y, xtype_col, j_col, k_col, cgmy_comp);
                }
                else {
                    tmp += _singular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, 0., xtype_col, j_col, k_col, true, cgmy_comp);
                    tmp += _singular_quadrature_dpsi_vs_CGMYkernel(x_ast, 0., b_y, xtype_col, j_col, k_col, true, cgmy_comp);
                }
            }
            else if ((a_y == 0.) && (b_y >0.)) {
                if (cgmy.Y<1) {
                    tmp += _nonsingular_quadrature_dpsi_vs_CGMYkernel(x_ast, 0., b_y, xtype_col, j_col, k_col, cgmy_comp);
                }
                else {
                    tmp += _singular_quadrature_dpsi_vs_CGMYkernel(x_ast, 0., b_y, xtype_col, j_col, k_col, true, cgmy_comp);
                }
            }
            else if ((a_y < 0.) && (b_y == 0.)) {
                if (cgmy.Y<1) {
                    tmp += _nonsingular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, 0, xtype_col, j_col, k_col, cgmy_comp);
                }
                else {
                    tmp += _singular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, 0., xtype_col, j_col, k_col, true, cgmy_comp);
                }
            }
            else {
                if (fabs(a_y)<0.01 || fabs(b_y)<0.01) {
                    tmp += _singular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, b_y, xtype_col, j_col, k_col, false, cgmy_comp);
                }
                else {
                    tmp += _nonsingular_quadrature_dpsi_vs_CGMYkernel(x_ast, a_y, b_y, xtype_col, j_col, k_col, cgmy_comp);
                }
            }
        }
        if (xtype_row == XBSpline) {
            ret += _weights(order,i) * tmp * d_phi(x_ast,j_row,k_row);
        }
        else {
            ret += _weights(order,i) * tmp * d_psi(x_ast,j_row,k_row);
        }
    }

    ret *= 0.5*(b-a);
    return ret;
}


template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::_nonsingular_quadrature_dpsi_vs_CGMYkernel(T x_ast, T a, T b, XType xtype_col,
                                                                    int j_col, int k_col,
                                                                    const CGMYUtils<T> &cgmy_comp) const
{
    T ret = 0.;
    if (fabs(a-b)<1e-16) return 0.;

    for (int i=1; i<=order; ++i) {
        T y = 0.5*(b-a)*_knots(order,i)+0.5*(b+a);
        if (xtype_col == XBSpline) {
            ret += _weights(order,i) * cgmy_comp.SecondTailIntegral(y) * d_phi(x_ast+y,j_col,k_col);
        }
        else {
            ret += _weights(order,i) * cgmy_comp.SecondTailIntegral(y) * d_psi(x_ast+y,j_col,k_col); ;
        }
    }
    return 0.5*(b-a)*ret;

}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::_singular_quadrature_dpsi_vs_CGMYkernel(T x_ast, T a, T b, XType xtype_col,
                                                                 int j_col, int k_col, bool cutoff,
                                                                 const CGMYUtils<T> &cgmy_comp) const
{
    T ret = 0.;
    assert((a>=0. && b>=0.) || (a<=0. && b<=0.));

    int q_j = 1;

    for (int j=0; j<=n+1; ++j) {
        T left, right;
        if (a>=0.) {
            if (j==0)   {
                if (cutoff) continue;
                left  = a;
                right = a+(b-a)*std::pow(sigma,n);
            }
            else if (j==n+1) continue;
            else        {
                left  = a+(b-a)*std::pow(sigma,n+1-j);
                right = a+(b-a)*std::pow(sigma,n-j);
            }
            q_j = std::max(1,int(std::ceil(mu*j)));
        }
        else {
            if (j==n+1)   {
                if (cutoff) continue;
                right = b;
                left  = b-(b-a)*std::pow(sigma,n);
            }
            else if (j==0) continue;
            else {
                right = b-(b-a)*std::pow(sigma,j);
                left  = b-(b-a)*std::pow(sigma,j-1);
            }
            q_j = std::max(1,int(std::ceil(mu*(n+1-j))));
        }
        T tmp = 0.;
        for (int i=1; i<=q_j; ++i) {
            T y = 0.5*(right-left)*_knots(q_j,i)+0.5*(right+left);
            if (xtype_col == XBSpline) {
                tmp += _weights(q_j,i) * cgmy_comp.SecondTailIntegral(y) * d_phi(x_ast+y,j_col,k_col);
            }
            else {
                tmp += _weights(q_j,i) * cgmy_comp.SecondTailIntegral(y) * d_psi(x_ast+y,j_col,k_col); ;
            }
        }
        ret += 0.5*(right-left)*tmp;
    }
    return ret;
}

}    //namespace lawa
