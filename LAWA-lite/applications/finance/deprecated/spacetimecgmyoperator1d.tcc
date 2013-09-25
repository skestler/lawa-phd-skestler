namespace lawa {

template <typename T, typename Basis, typename CGMYOperator>
SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::SpaceTimeCGMYOperator1D(const Basis& _basis, const CGMYOperator &_cgmy_x)
    : basis(_basis),
      d_t(basis.first), id_t(basis.first), id_x(basis.second), cgmy_x(_cgmy_x),
      phi_t(basis.first.mra), d_phi_t(basis.first.mra, 1), phi_x(basis.second.mra),
      psi_t(basis.first), d_psi_t(basis.first, 1), psi_x(basis.second),
      integral_sfsf_t(phi_t, phi_t), d_integral_sfsf_t(phi_t, d_phi_t), integral_sfsf_x(phi_x, phi_x),
      integral_sfw_t(phi_t, psi_t),  d_integral_sfw_t(phi_t, d_psi_t),  integral_sfw_x(phi_x, psi_x),
      integral_wsf_t(psi_t, phi_t),  d_integral_wsf_t(psi_t, d_phi_t),  integral_wsf_x(psi_x, phi_x),
      integral_ww_t(psi_t, psi_t),   d_integral_ww_t(psi_t, d_psi_t),   integral_ww_x(psi_x, psi_x)
{
}

template <typename T, typename Basis, typename CGMYOperator>
T
SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::getc() const
{
    return 1.;
}

template <typename T, typename Basis, typename CGMYOperator>
const Basis&
SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis, typename CGMYOperator>
T
SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::operator()(XType row_xtype_t, int j_row_t, int k_row_t,
                                              XType row_xtype_x, int j_row_x, int k_row_x,
                                              XType col_xtype_t,int j_col_t, int k_col_t,
                                              XType col_xtype_x,int j_col_x, int k_col_x) const
{
    T val_t = 0;
    T d_val_t = 0;
    T val_x = 0;
    T cgmy_val_x = 0;

    if(row_xtype_t == XBSpline){
         if(col_xtype_t == XBSpline){
             val_t = integral_sfsf_t(j_row_t, k_row_t, j_col_t, k_col_t);
             d_val_t = d_integral_sfsf_t(j_row_t, k_row_t, j_col_t, k_col_t);
         }
         else{
             val_t = integral_sfw_t(j_row_t, k_row_t, j_col_t, k_col_t);
             d_val_t = d_integral_sfw_t(j_row_t, k_row_t, j_col_t, k_col_t);
         }
    }
    else{
         if(col_xtype_t == XBSpline){
             val_t = integral_wsf_t(j_row_t, k_row_t, j_col_t, k_col_t);
             d_val_t = d_integral_wsf_t(j_row_t, k_row_t, j_col_t, k_col_t);
         }
         else{
             val_t = integral_ww_t(j_row_t, k_row_t, j_col_t, k_col_t);
             d_val_t = d_integral_ww_t(j_row_t, k_row_t, j_col_t, k_col_t);
         }
    }

    if(row_xtype_x == XBSpline){
         if(col_xtype_x == XBSpline){
             val_x = integral_sfsf_x(j_row_x, k_row_x, j_col_x, k_col_x);
         }
         else{
             val_x = integral_sfw_x(j_row_x, k_row_x, j_col_x, k_col_x);
         }
    }
    else{
         if(col_xtype_x == XBSpline){
             val_x = integral_wsf_x(j_row_x, k_row_x, j_col_x, k_col_x);
         }
         else{
             val_x = integral_ww_x(j_row_x, k_row_x, j_col_x, k_col_x);
         }
    }

    cgmy_val_x = cgmy_x(row_xtype_x, j_row_x, k_row_x, col_xtype_x, j_col_x, k_col_x);

    return d_val_t * val_x  + val_t * cgmy_val_x;

}

template <typename T, typename Basis, typename CGMYOperator>
T
SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::operator()(const Index2D &row_index,
                                              const Index2D &col_index) const
{
    return operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                      row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                      col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                      col_index.index2.xtype, col_index.index2.j, col_index.index2.k);

}

}    //namespace lawa
