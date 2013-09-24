namespace lawa {

template <typename T, typename Basis2D>
OptimizedH1Preconditioner2D<T,Basis2D>::OptimizedH1Preconditioner2D(const Basis2D &basis,
                                                                    T _a_x, T _a_y, T _c)
    : basis_x(basis.first), basis_y(basis.second), a_x(_a_x), a_y(_a_y), c(_c),
      factor_a_x(_a_x), factor_a_y(_a_y), factor_c(_c)
{

}

template <typename T, typename Basis2D>
void
OptimizedH1Preconditioner2D<T,Basis2D>::setParameters(T _a_x, T _a_y, T _c)
{
    a_x = _a_x;
    a_y = _a_y;
    c   = _c;
    factor_a_x = _a_x;
    factor_a_y = _a_y;
    factor_c = _c;
}

template <typename T, typename Basis2D>
void
OptimizedH1Preconditioner2D<T,Basis2D>::setThetaTimeStepParameters(T theta, T timestep)
{
    factor_a_x = a_x*theta*timestep;
    factor_a_y = a_y*theta*timestep;
}

template <typename T, typename Basis2D>
T
OptimizedH1Preconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, long k1,
                                                   XType xtype2, int j2, long k2) const
{
    T id_x=0., dd_x=0., id_y=0., dd_y=0.;
    id_x = basis_x.generator(xtype1).getL2Norm(j1,k1);
    id_x *= id_x;
    dd_x = basis_x.generator(xtype1).getH1SemiNorm(j1,k1);
    dd_x *= dd_x;
    id_y = basis_y.generator(xtype2).getL2Norm(j2,k2);
    id_y *= id_y;
    dd_y = basis_y.generator(xtype2).getH1SemiNorm(j2,k2);
    dd_y *= dd_y;

    return (T)1./std::sqrt(factor_a_x*dd_x*id_y + factor_a_y*id_x*dd_y + factor_c*id_x*id_y);
}

template <typename T, typename Basis2D>
T
OptimizedH1Preconditioner2D<T,Basis2D>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k);
}

template <typename T, typename Basis2D>
T
OptimizedH1Preconditioner2D<T,Basis2D>::operator[](const Index2D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k);
}

}   // namespace lawa
