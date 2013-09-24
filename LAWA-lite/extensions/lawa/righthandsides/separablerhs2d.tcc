namespace lawa {

template<typename T, typename Basis2D>
SeparableRHS2D<T, Basis2D>::SeparableRHS2D
                            (const Basis2D& _basis, const SeparableFunction2D<T>& _F,
                             const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_x,
                             const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_y,
                             int order)
    : basis(_basis), F(_F), deltas_x(_deltas_x), deltas_y(_deltas_y),
      integralf_x(F.F_x, basis.first), integralf_y(F.F_y, basis.second)
{
    integralf_x.quadrature.setOrder(order);
    integralf_y.quadrature.setOrder(order);
}

template<typename T, typename Basis2D>
T
SeparableRHS2D<T, Basis2D>::operator()(XType xtype_x, int j_x, long k_x,
                                       XType xtype_y, int j_y, long k_y) const
{
    T val_x = 0;
    T val_y = 0;

    val_x += integralf_x(j_x, k_x, xtype_x, 0);
    for (int i=1; i<=deltas_x.numRows(); ++i) {
        val_x += deltas_x(i,2) * basis.first.generator(xtype_x)(deltas_x(i,1),j_x,k_x,0);
    }

    val_y += integralf_y(j_y, k_y, xtype_y, 0);
    for (int i=1; i<=deltas_y.numRows(); ++i) {
        val_y += deltas_y(i,2) * basis.second.generator(xtype_y)(deltas_y(i,1),j_y,k_y,0);
    }

    return val_x * val_y;
}

template<typename T, typename Basis2D>
T
SeparableRHS2D<T, Basis2D>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

} // namespace lawa

