namespace lawa {

template <typename T, typename Basis1D>
FinanceOperator1D<T, CGMY, Basis1D>::FinanceOperator1D
                                    (const Basis1D& _basis,
                                     const ProcessParameters1D<T,CGMY> &_processparameters,
                                     T _R1, T _R2, int order, const int _internal_compression_level,
                                     T _convection, T _reaction, const bool _use_predef_convection,
                                     const bool _use_predef_reaction)
    : basis(_basis), processparameters(_processparameters), kernel(processparameters),
      R1(_R1), R2(_R2), internal_compression_level(_internal_compression_level),
      convection(_convection), reaction(_reaction),
      use_predef_convection(_use_predef_convection),
      use_predef_reaction(_use_predef_reaction),
      OneDivSqrtR2pR1(1./std::sqrt(R2+R1)), OneDivR2pR1(1./(R2+R1)), R1DivR1pR2(R1/(R1+R2)),
      integral(basis,basis)
{
    if (basis.d==3) {
        bool is_realline_basis=flens::IsSame<Basis1D, Basis<T,Primal,R,CDF> >::value;
        assert(is_realline_basis);
    }
    if (use_predef_convection) {
        convection = kernel.ExpXmOnemX_k;
        reaction   = 0.;
    }

}

template <typename T, typename Basis1D>
T
FinanceOperator1D<T, CGMY, Basis1D>::operator()(XType xtype1, int j1, int k1,
                                                XType xtype2, int j2, int k2) const
{
   typedef typename std::map<T,T>::const_iterator const_it;

   if (internal_compression_level>-1) {
       //T compr_c=0.1;
       T compr_c=1.;
       T compr_alpha=T(2*basis.d)/T(2*basis.d_+processparameters.k_Y);
       if (xtype1==XWavelet && xtype2==XWavelet) {
           int J=internal_compression_level;
           T delta = compr_c*std::max(pow2i<T>(-std::min(j1,j2)),
                                      pow2i<T>(-(J-1)+compr_alpha*(2*(J-1)-j1-j2)));
           if ( ( distance(basis.generator(xtype1).support(j1,k1),
                           basis.generator(xtype2).support(j2,k2)) > delta/(R1+R2))
                && (k1 > basis.rangeJL(j1).lastIndex())
                && (k1 < basis.rangeJR(j1).firstIndex())
                && (k2 > basis.rangeJL(j2).lastIndex())
                && (k2 < basis.rangeJR(j2).firstIndex()) ) return 0.;
       }
   }

   T pde_val = 0.;
   pde_val += convection * OneDivR2pR1 * integral(j1,k1,xtype1,0,j2,k2,xtype2,1);
   if (!use_predef_reaction) {
       pde_val += reaction * integral(j1,k1,xtype1,0,j2,k2,xtype2,0);
   }


   T int_val = 0.;

   GeMatrix<FullStorage<T,ColMajor> > varphi_row_deltas, varphi_col_deltas;
   varphi_row_deltas = computeDeltas<T,Basis1D>(basis,j1,k1,xtype1);
   varphi_col_deltas = computeDeltas<T,Basis1D>(basis,j2,k2,xtype2);

   if (R1!=0 && R2!=1) {
       varphi_row_deltas(_,1)  *=(R1+R2);  varphi_row_deltas(_,1)-=R1;
       varphi_row_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
       varphi_col_deltas(_,1)  *=(R1+R2);  varphi_col_deltas(_,1)-=R1;
       varphi_col_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
   }

   T part1=0., part2=0.;
   if (basis.d==2) {
       for (int lambda=varphi_row_deltas.rows().firstIndex();
                lambda<=varphi_row_deltas.rows().lastIndex(); ++lambda) {

           T x = varphi_row_deltas(lambda,1);
           if (fabs(varphi_row_deltas(lambda,2)) < 1e-14) continue;

           part1 += OneDivSqrtR2pR1*varphi_row_deltas(lambda,2)
                    *basis.generator(xtype2)((x+R1)* OneDivR2pR1,j2,k2,0)*kernel.c3;


           for (int mu=varphi_col_deltas.rows().firstIndex();
                    mu<=varphi_col_deltas.rows().lastIndex(); ++mu) {
               if (fabs(varphi_col_deltas(mu,2)) < 1e-14) continue;

               T y = varphi_col_deltas(mu,1);
               T c = varphi_col_deltas(mu,2)*varphi_row_deltas(lambda,2);


               if (fabs(x-y)>1e-10)  {
                   T value_tailintegral;
                   const_it it   = values_tailintegral.find(y-x);
                   const_it last = values_tailintegral.end();
                   if (it != last) value_tailintegral = (*it).second;
                   else {
                       value_tailintegral=kernel.ForthTailIntegral(y-x);
                       values_tailintegral[y-x] = value_tailintegral;
                   }
                   if (y-x>0)  part2 += c * (value_tailintegral - kernel.constants[2]);
                   else        part2 += c * (value_tailintegral - kernel.constants[3]);
               }
           }
       }
       int_val = part1 + part2;
   }

   else if (basis.d==3) {

       int_val -= kernel.c3*integral(j1,k1,xtype1,1,j2,k2,xtype2,1);

       for (int mu=varphi_col_deltas.rows().firstIndex();
                mu<=varphi_col_deltas.rows().lastIndex(); ++mu) {

           T y = varphi_col_deltas(mu,1);
           if (fabs(varphi_col_deltas(mu,2)) < 1e-14) continue;

           part1 += varphi_col_deltas(mu,2)*basis.generator(xtype1)(y,j1,k1,0)*kernel.c4;
           part1 -= varphi_col_deltas(mu,2)*basis.generator(xtype1)(y,j1,k1,1)*kernel.c5;

           for (int lambda=varphi_row_deltas.rows().firstIndex();
                    lambda<=varphi_row_deltas.rows().lastIndex(); ++lambda) {

               if (fabs(varphi_row_deltas(lambda,2)) < 1e-14) continue;
               T x = varphi_row_deltas(lambda,1);
               T c = varphi_col_deltas(mu,2)*varphi_row_deltas(lambda,2);
               if (x!=y)  {
                   if (y-x>0)  {
                       part2 -= c * (kernel.SixthTailIntegral(y-x) - kernel.constants[6]);
                   }
                   else    {
                       part2 -= c * (kernel.SixthTailIntegral(y-x) - kernel.constants[7]);
                   }
               }
           }
       }
       int_val += part1 + part2;
   }
   else {
       assert(0);
   }
   return pde_val - int_val;
}

}    //namespace lawa
