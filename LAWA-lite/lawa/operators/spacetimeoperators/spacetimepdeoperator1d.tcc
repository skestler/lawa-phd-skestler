namespace lawa{
    
template <typename T, typename Basis>    
SpaceTimePDEOperator1D<T, Basis>::SpaceTimePDEOperator1D(const Basis& _basis, const T _diffusion, 
                                                         const T _convection, const T _reaction,
                                                         const T _timederivfactor)
    : basis(_basis), diffusion(_diffusion), convection(_convection), reaction(_reaction),
      timederivfactor(_timederivfactor),
      integral_t(_basis.first, _basis.first), integral_x(_basis.second, _basis.second)
{
}
    
template <typename T, typename Basis>      
T
SpaceTimePDEOperator1D<T, Basis>::operator()(XType row_xtype_t, int j1_t, long k1_t,
                                             XType row_xtype_x, int j1_x, long k1_x,
                                             XType col_xtype_t, int j2_t, long k2_t,
                                             XType col_xtype_x, int j2_x, long k2_x) const
{
    T val_t =    integral_t(j1_t, k1_t, row_xtype_t, 0, j2_t, k2_t, col_xtype_t, 0);
    T d_val_t =  integral_t(j1_t, k1_t, row_xtype_t, 0, j2_t, k2_t, col_xtype_t, 1);
    T val_x =    integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0);
    T d_val_x =  integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 1);
    T dd_val_x = integral_x(j1_x, k1_x, row_xtype_x, 1, j2_x, k2_x, col_xtype_x, 1);
    
    // (v1 * u1_t)*(v2*u2) + diffusion * (v1 * u1)*(v2_x*u2_x) 
    //  + convection * (v1 * u1)*(v2*u2_x) + reaction * (v1 * u1)*(v2 * u2)
    return timederivfactor * d_val_t * val_x  + diffusion * val_t * dd_val_x 
            + convection * val_t * d_val_x + reaction * val_t * val_x;

}

template <typename T, typename Basis>  
T
SpaceTimePDEOperator1D<T, Basis>::operator()(XType xtype_t, int j_t, long k_t,
                                             XType xtype_x, int j_x, long k_x) const
{
    return operator()(xtype_t, j_t, k_t, xtype_x, j_x, k_x, 
                      xtype_t, j_t, k_t, xtype_x, j_x, k_x);
}

template <typename T, typename Basis>  
T
SpaceTimePDEOperator1D<T, Basis>::operator()(const Index2D &row_index, 
                                              const Index2D &col_index) const
{
    return operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                      row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                      col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                      col_index.index2.xtype, col_index.index2.j, col_index.index2.k);

}
    
} // namespace lawa
