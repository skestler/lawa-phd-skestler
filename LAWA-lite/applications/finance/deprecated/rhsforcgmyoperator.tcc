namespace lawa {

template <typename T, typename Basis>
RHSForCGMYOperator1D<T,Basis>::RHSForCGMYOperator1D(const Basis &_basis, T (*_f)(T), T (*_df)(T),
                                                const DenseVector<Array<T> > &_singularPoints,
                                                T C, T G, T M, T Y, int _order)
    : basis(_basis),  phi(basis.mra), psi(basis),
      f(_f,_singularPoints), df(_df,_singularPoints), cgmy(C,G,M,Y), order(_order)
{
    int jmin = basis.j0;
    int d    = basis.d;
    T step = 1.0/(1<<(jmin+7));
    Support<T> supp;
    phi_deltas.engine().resize(phi.singularSupport(jmin,0).length(),2);
    phi_deltas(_,1) = phi.singularSupport(jmin,0);
    supp = phi.support(jmin,0);

    PrimalSpline  dM1_th_deriv_phi(basis.mra,d-1);
    PrimalWavelet dM1_th_deriv_psi(basis,d-1);

    for (int i = 1; i<=phi_deltas.numRows(); ++i) {
        phi_deltas(i,2) = ((phi_deltas(i,1)==supp.l2) ? 0.0 : dM1_th_deriv_phi(std::min(phi_deltas(i,1)+step, supp.l2),jmin,0))
                        - ((phi_deltas(i,1)==supp.l1) ? 0.0 : dM1_th_deriv_phi(std::max(phi_deltas(i,1)-step, supp.l1),jmin,0));

    }

    psi_deltas.engine().resize(psi.singularSupport(jmin,0).length(),2);
    psi_deltas(_,1) = psi.singularSupport(jmin,0);
    supp = psi.support(jmin,0);

    for (int i = 1; i<=psi_deltas.numRows(); ++i) {
        psi_deltas(i,2) = ((psi_deltas(i,1)==supp.l2) ? 0.0 : dM1_th_deriv_psi(std::min(psi_deltas(i,1)+step, supp.l2),jmin,0))
                        - ((psi_deltas(i,1)==supp.l1) ? 0.0 : dM1_th_deriv_psi(std::max(psi_deltas(i,1)-step, supp.l1),jmin,0));
    }


    T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(_order, _order);
    _weights.engine().resize(_order, _order);
    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=order; ++k) {
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
RHSForCGMYOperator1D<T,Basis>::operator()(XType xtype, int j, int k) const
{
    T ret = 0.;
    GeMatrix<FullStorage<T,ColMajor> > rhs_psi_deltas=RHSForCGMYOperator1D<T,Basis>::_computeDeltas(xtype,j,k);
    if (basis.d == 2) {
        for (int i=rhs_psi_deltas.rows().firstIndex(); i<=rhs_psi_deltas.rows().lastIndex(); ++i) {
            T x = rhs_psi_deltas(i,1), v = rhs_psi_deltas(i,2);
            T tmp = RHSForCGMYOperator1D<T,Basis>::_evaluate_acal_u(x,-80,80);
            ret += v * tmp;
            //std::cout << Index1D(j,k,xtype) << ": " << x << " " << v << ": " << tmp << std::endl;
        }
    }
    else {
        assert(0);
    }
    return ret;
}

template <typename T, typename Basis>
T
RHSForCGMYOperator1D<T,Basis>::operator()(const Index1D &lambda) const
{
    return RHSForCGMYOperator1D<T,Basis>::operator()(lambda.xtype, lambda.j, lambda.k);
}

template <typename T, typename Basis>
GeMatrix<FullStorage<T,ColMajor> >
RHSForCGMYOperator1D<T,Basis>::_computeDeltas(XType xtype, int j, int k) const
{
    int jmin=basis.j0;
    int d   =basis.d;
    GeMatrix<FullStorage<T,ColMajor> > ret;
    if (xtype == XBSpline) {
        ret.engine().resize(phi_deltas.numRows(),phi_deltas.numCols());
        for (int i = 1; i<=ret.numRows(); ++i) {
            ret(i,1) = phi_deltas(i,1)+pow2i<T>(-jmin)*k;
            ret(i,2) = phi_deltas(i,2);
        }
    }
    else {
        ret.engine().resize(psi_deltas.numRows(),psi_deltas.numCols());
        for (int i = 1; i<=ret.numRows(); ++i) {
            ret(i,1) = pow2i<T>(jmin-j)*psi_deltas(i,1)+pow2i<T>(-j)*k;
            ret(i,2) = pow2i<T>((d-1)*(j-jmin))*pow2ih<T>(j-jmin)*psi_deltas(i,2);
        }
    }
    return ret;
}

template <typename T, typename Basis>
T
RHSForCGMYOperator1D<T,Basis>::_evaluate_acal_u(T x_ast, T a, T b) const
{
    flens::DenseVector<Array<T> > singularPoints_integration(3);
    singularPoints_integration = a, 0., b;
    flens::DenseVector<Array<T> > singularPoints_function(f.singularPoints.length());
    for (int i=1; i<=singularPoints_function.length(); ++i) {
        singularPoints_function(i) = f.singularPoints(i) - x_ast;
    }

    DenseVector<Array<T> > tempPoints(3 + df.singularPoints.length());
    std::merge(singularPoints_integration.engine().data(),
               singularPoints_integration.engine().data() + 3,
               singularPoints_function.engine().data(),
               singularPoints_function.engine().data() + singularPoints_function.length(),
               tempPoints.engine().data());


    T ret = 0.0;
    ret += cgmy.c3 * f(x_ast);
    for (int k1=tempPoints.firstIndex(); k1<tempPoints.lastIndex(); ++k1) {
        T h = tempPoints(k1+1)-tempPoints(k1);
        T max_interval_length = 1.;
        //if ((tempPoints(k1)>f.singularPoints(1)-3) && (tempPoints(k1)<f.singularPoints(1)+3)) {
        //    max_interval_length = 1;
        //}
        int NumOfSubInt = 1;
        do {
            h *= 0.5; NumOfSubInt *= 2;
        } while(h > max_interval_length);
        T x = tempPoints(k1);
        for (int i=0; i<NumOfSubInt; ++i) {
            T tmp = _quadrature_du_vs_CGMYkernel(x_ast, x+i*h, x+(i+1)*h);
            //std::cout << "x = " << x_ast << ", [" << x+i*h << ", " << x+(i+1)*h << "]:" << tmp << std::endl;
            ret += tmp;
        }
    }
    return ret;
}

template <typename T, typename Basis>
T
RHSForCGMYOperator1D<T,Basis>::_quadrature_du_vs_CGMYkernel(T x_ast, T a, T b) const
{
    T ret = 0.;
    T c1 = (b-a)/2., d1 = (b+a)/2.;
    for (int i = 1; i <= order; ++i) {
        T x_tmp = c1*_knots(order,i)+d1;
        //std::cout << "   " << x_tmp << " " <<  _weights(order,i) << " " << df(x_ast+x_tmp) << " " << cgmy.ThirdTailIntegral(x_tmp) << std::endl;
        ret += _weights(order,i)*df(x_ast+x_tmp)*cgmy.ThirdTailIntegral(x_tmp);
    }
    return 0.5*fabs((b-a))*ret;
}

}    //namspace lawa
