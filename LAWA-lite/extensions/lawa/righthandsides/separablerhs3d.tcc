namespace lawa {

template<typename T, typename Basis3D>
SeparableRHS3D<T, Basis3D>::SeparableRHS3D
                            (const Basis3D& _basis, const SeparableFunction3D<T>& _F,
                             const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_x,
                             const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_y,
                             const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_z,
                             int order)
    : basis(_basis), F(_F),
      deltas_x(_deltas_x), deltas_y(_deltas_y), deltas_z(_deltas_z),
      integralf_x(F.F_x, basis.first), integralf_y(F.F_y, basis.second),
      integralf_z(F.F_z, basis.third)
{
    integralf_x.quadrature.setOrder(order);
    integralf_y.quadrature.setOrder(order);
    integralf_z.quadrature.setOrder(order);
}

template<typename T, typename Basis3D>
T
SeparableRHS3D<T, Basis3D>::operator()(XType xtype_x, int j_x, long k_x,
                                       XType xtype_y, int j_y, long k_y,
                                       XType xtype_z, int j_z, long k_z) const
{
    T val_x = 0;
    T val_y = 0;
    T val_z = 0;

    val_x = integralf_x(j_x, k_x, xtype_x, 0);
    for (int i=1; i<=deltas_x.numRows(); ++i) {
        val_x += deltas_x(i,2) * basis.first.generator(xtype_x)(deltas_x(i,1),j_x,k_x,0);
    }

    val_y = integralf_y(j_y, k_y, xtype_y, 0);
    for (int i=1; i<=deltas_y.numRows(); ++i) {
        val_y += deltas_y(i,2) * basis.second.generator(xtype_y)(deltas_y(i,1),j_y,k_y,0);
    }

    val_z = integralf_z(j_z, k_z, xtype_z, 0);
    for (int i=1; i<=deltas_z.numRows(); ++i) {
        val_z += deltas_z(i,2) * basis.third.generator(xtype_z)(deltas_z(i,1),j_z,k_z,0);
    }

    return val_x * val_y * val_z;
}

template<typename T, typename Basis3D>
T
SeparableRHS3D<T, Basis3D>::operator()(const Index3D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k,
                            index.index3.xtype, index.index3.j, index.index3.k);
}

} // namespace lawa

