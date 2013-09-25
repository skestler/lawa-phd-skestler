namespace lawa {

template <typename T, typename Basis>
CGMYOperator1D_Fast<T, Basis>::CGMYOperator1D_Fast(const Basis& _basis, T _diffusion, T _convection,
                                                   T _reaction, T _C, T _G, T _M, T _Y)
    : basis(_basis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      C(_C), G(_G), M(_M), Y(_Y),
      cgmy(C,G,M,Y),
      phi(basis.mra), d_phi(basis.mra, 1), delta_phi(basis.mra,basis.d-1),
      psi(basis), d_psi(basis, 1), delta_psi(basis,basis.d-1),
      integral_sfsf(phi, phi), d_integral_sfsf(phi, d_phi), dd_integral_sfsf(d_phi, d_phi),
      integral_sfw(phi, psi),  d_integral_sfw(phi, d_psi),  dd_integral_sfw(d_phi, d_psi),
      integral_wsf(psi, phi),  d_integral_wsf(psi, d_phi),  dd_integral_wsf(d_psi, d_phi),
      integral_ww(psi,psi),    d_integral_ww(psi, d_psi),   dd_integral_ww(d_psi,d_psi)
{
    assert(diffusion>=0.);
    assert(reaction >=0.);
    assert(C>=0.);

}

template <typename T, typename Basis>
T
CGMYOperator1D_Fast<T,Basis>::getc() const
{
    return reaction;
}

template <typename T, typename Basis>
const Basis&
CGMYOperator1D_Fast<T,Basis>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis>
GeMatrix<FullStorage<T,ColMajor> >
CGMYOperator1D_Fast<T, Basis>::computeDeltas(XType xtype, int j, int k) const
{

    GeMatrix<FullStorage<T,ColMajor> > ret;

    if (xtype == XBSpline) {
        ret.engine().resize(phi.singularSupport(j,k).length(),2);
        ret(_,1) = phi.singularSupport(j,k);
        Support<T> supp = phi.support(j,k);
        T step = pow2i<T>(-(j+5)); // 1.0/(1<<(j+1));
        for (int i = 1; i<=ret.numRows(); ++i) {
            ret(i,2) = ((ret(i,1)==supp.l2) ? 0.0 : delta_phi(std::min(ret(i,1)+step, supp.l2),j,k))
                     - ((ret(i,1)==supp.l1) ? 0.0 : delta_phi(std::max(ret(i,1)-step, supp.l1),j,k));
        }
    }
    else {
        ret.engine().resize(psi.singularSupport(j,k).length(),2);
        ret(_,1) = psi.singularSupport(j,k);
        Support<T> supp = psi.support(j,k);
        T step = pow2i<T>(-(j+5)); // 1.0/(1<<(j+1));
        for (int i = 1; i<=ret.numRows(); ++i) {
            ret(i,2) = ((ret(i,1)==supp.l2) ? 0.0 : delta_psi(std::min(ret(i,1)+step, supp.l2),j,k))
                     - ((ret(i,1)==supp.l1) ? 0.0 : delta_psi(std::max(ret(i,1)-step, supp.l1),j,k));
        }
    }
    return ret;
}

template <typename T, typename Basis>
T
CGMYOperator1D_Fast<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                          XType xtype2, int j2, int k2) const
{
    //PDE part
    T val = 0;
    T d_val = 0;
    T dd_val = 0;

    if(xtype1 == XBSpline){
         if(xtype2 == XBSpline){
             val = integral_sfsf(j1, k1, j2, k2);
             d_val = d_integral_sfsf(j1, k1, j2, k2);
             dd_val = dd_integral_sfsf(j1, k1, j2, k2);
         }
         else{
             val = integral_sfw(j1, k1, j2, k2);
             d_val = d_integral_sfw(j1, k1, j2, k2);
             dd_val = dd_integral_sfw(j1, k1, j2, k2);
         }
    }
    else{
         if(xtype2 == XBSpline){
             val = integral_wsf(j1, k1, j2, k2);
             d_val = d_integral_wsf(j1, k1, j2, k2);
             dd_val = dd_integral_wsf(j1, k1, j2, k2);
         }
         else{
             val = integral_ww(j1, k1, j2, k2);
             d_val = d_integral_ww(j1, k1, j2, k2);
             dd_val = dd_integral_ww(j1, k1, j2, k2);
         }
    }

    GeMatrix<FullStorage<T,ColMajor> > psi_row_deltas = computeDeltas(xtype1,j1,k1);
    GeMatrix<FullStorage<T,ColMajor> > psi_col_deltas = computeDeltas(xtype2,j2,k2);

    T int_val = 0.;

    T part1=0., part2=0.;
    if (basis.d==2) {

        for (int lambda=psi_row_deltas.rows().firstIndex(); lambda<=psi_row_deltas.rows().lastIndex(); ++lambda) {
            T x = psi_row_deltas(lambda,1);
            if (fabs(psi_row_deltas(lambda,2)) < 1e-14) continue;


            if (xtype2 == XBSpline)    part1 += psi_row_deltas(lambda,2)*phi(x,j2,k2)*cgmy.c3;
            else                    part1 += psi_row_deltas(lambda,2)*psi(x,j2,k2)*cgmy.c3;

            for (int mu=psi_col_deltas.rows().firstIndex(); mu<=psi_col_deltas.rows().lastIndex(); ++mu) {
                if (fabs(psi_col_deltas(mu,2)) < 1e-14) continue;
                T y = psi_col_deltas(mu,1);
                T c = psi_col_deltas(mu,2)*psi_row_deltas(lambda,2);
                if (fabs(x-y)>1e-10)  {
                    if (y-x>0)    part2 += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[2]);
                    else         part2 += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[3]);
                }
            }
        }
        int_val = part1 + part2;
    }

    else if (basis.d==3) {

        if         ((xtype1 == XBSpline) && (xtype2 == XBSpline))     int_val -= cgmy.c3*dd_integral_sfsf(j1,k1,j2,k2);
        else if ((xtype1 == XBSpline) && (xtype2 == XWavelet))    int_val -= cgmy.c3*dd_integral_sfw(j1,k1,j2,k2);
        else if ((xtype1 == XWavelet) && (xtype2 == XBSpline))    int_val -= cgmy.c3*dd_integral_wsf(j1,k1,j2,k2);
        else                                                    int_val -= cgmy.c3*dd_integral_ww(j1,k1,j2,k2);

        for (int mu=psi_col_deltas.rows().firstIndex(); mu<=psi_col_deltas.rows().lastIndex(); ++mu) {
            T y = psi_col_deltas(mu,1);
            if (fabs(psi_col_deltas(mu,2)) < 1e-14) continue;

            if (xtype1 == XBSpline)    {
                part1 += psi_col_deltas(mu,2)*phi(y,j1,k1)*cgmy.c4;
                part1 -= psi_col_deltas(mu,2)*d_phi(y,j1,k1)* cgmy.c5;
            }
            else {
                part1 += psi_col_deltas(mu,2)*psi(y,j1,k1)*cgmy.c4;
                part1 -= psi_col_deltas(mu,2)*d_psi(y,j1,k1)* cgmy.c5;
            }

            for (int lambda=psi_row_deltas.rows().firstIndex(); lambda<=psi_row_deltas.rows().lastIndex(); ++lambda) {
                if (fabs(psi_row_deltas(lambda,2)) < 1e-14) continue;
                T x = psi_row_deltas(lambda,1);
                T c = psi_col_deltas(mu,2)*psi_row_deltas(lambda,2);
                if (x!=y)  {
                    if (y-x>0)    {
                        part2 -= c * (cgmy.SixthTailIntegral(y-x) - cgmy.constants[6]);
                    }
                    else     {
                        part2 -= c * (cgmy.SixthTailIntegral(y-x) - cgmy.constants[7]);
                    }
                }
            }
        }
        int_val += part1 + part2;
    }
    else {
        assert(0);
    }
    return diffusion*dd_val + convection * d_val +  reaction * val - int_val;
}

template <typename T, typename Basis>
T
CGMYOperator1D_Fast<T, Basis>::operator()(XType xtype_row, int j_row, int k_row,
                                          XType xtype_col, int j_col, int k_col, T R1, T R2) const
{
    //only to be used for interval discretizations
    assert(diffusion==0);
    assert(reaction==0);

    T OneDivSqrtR2pR1 = 1./sqrt(R2+R1);
    T OneDivR2pR1 = 1./(R2+R1);

    GeMatrix<FullStorage<T,ColMajor> > psi_row_deltas = computeDeltas(xtype_row,j_row,k_row);
    GeMatrix<FullStorage<T,ColMajor> > psi_col_deltas = computeDeltas(xtype_col,j_col,k_col);

    psi_row_deltas(_,1)  *=(R1+R2);  psi_row_deltas(_,1)-=R1;
    psi_row_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
    psi_col_deltas(_,1)  *=(R1+R2);  psi_col_deltas(_,1)-=R1;
    psi_col_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;

    T int_val = 0.;

    //std::cout << "(" << j_row << ", " << k_row << "): " << psi_row_deltas << std::endl;
    //std::cout << "(" << j_col << ", " << k_col << "): " << psi_col_deltas << std::endl;

    if (basis.d==2) {
        T part1 = 0.;
        T part2 = 0.;

        for (int lambda=psi_row_deltas.rows().firstIndex(); lambda<=psi_row_deltas.rows().lastIndex(); ++lambda) {
            T x = psi_row_deltas(lambda,1);

            if (fabs(psi_row_deltas(lambda,2)) < 1e-14) continue;

            if (xtype_col == XBSpline)    part1 += OneDivSqrtR2pR1*psi_row_deltas(lambda,2)*phi((x+R1)* OneDivR2pR1,j_col,k_col)*cgmy.c3;
            else                        part1 += OneDivSqrtR2pR1*psi_row_deltas(lambda,2)*psi((x+R1)* OneDivR2pR1,j_col,k_col)*cgmy.c3;

            for (int mu=psi_col_deltas.rows().firstIndex(); mu<=psi_col_deltas.rows().lastIndex(); ++mu) {
                if (fabs(psi_col_deltas(mu,2)) < 1e-14) continue;
                T y = psi_col_deltas(mu,1);
                T c = psi_col_deltas(mu,2)*psi_row_deltas(lambda,2);
                if (fabs(x-y)>1e-10)  {
                    if (y-x>0)    {
                        part2 += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[2]);
                        //std::cout << y-x << " " << cgmy.ForthTailIntegral(y-x) << " " << cgmy.constants[2] << std::endl;
                    }
                    else {
                        part2 += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[3]);
                        //std::cout << y-x << " " << cgmy.ForthTailIntegral(y-x) << " " << cgmy.constants[3] << std::endl;
                    }
                }
            }
            //std::cout << "part1 = " << part1 << ", part2 = " << part2 << std::endl;
            //getchar();
        }
        int_val = part1 + part2;
    }

    else {
        std::cout << "Not stable for d>=3. Use CGMYOperator1D instead." << std::endl;
        assert(0);
    }
    return int_val;
}

template <typename T, typename Basis>
T
CGMYOperator1D_Fast<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return CGMYOperator1D_Fast<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                col_index.xtype, col_index.j, col_index.k);
}

}    //namespace lawa
