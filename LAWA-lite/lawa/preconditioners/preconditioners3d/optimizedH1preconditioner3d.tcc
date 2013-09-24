namespace lawa {

template <typename T, typename Basis3D>
OptimizedH1Preconditioner3D<T,Basis3D>::OptimizedH1Preconditioner3D(const Basis3D &basis,
                                                                    T _a_x, T _a_y, T _a_z, T _c)
    : basis_x(basis.first), basis_y(basis.second), basis_z(basis.third),
      a_x(_a_x), a_y(_a_y), a_z(_a_z), c(_c)
{

}

template <typename T, typename Basis3D>
void
OptimizedH1Preconditioner3D<T,Basis3D>::setParameters(T _a_x, T _a_y, T _a_z, T _c)
{
    a_x = _a_x; a_y = _a_y; a_z = _a_z; c = _c;
}

template <typename T, typename Basis3D>
T
OptimizedH1Preconditioner3D<T,Basis3D>::operator()(XType xtype_x, int j_x, long k_x,
                                                   XType xtype_y, int j_y, long k_y,
                                                   XType xtype_z, int j_z, long k_z) const
{
    T id_x=0., dd_x=0., id_y=0., dd_y=0., id_z=0., dd_z=0.;
    id_x = basis_x.generator(xtype_x).getL2Norm(j_x,k_x);
    id_x *= id_x;
    dd_x = basis_x.generator(xtype_x).getH1SemiNorm(j_x,k_x);
    dd_x *= dd_x;
    id_y = basis_y.generator(xtype_y).getL2Norm(j_y,k_y);
    id_y *= id_y;
    dd_y = basis_y.generator(xtype_y).getH1SemiNorm(j_y,k_y);
    dd_y *= dd_y;
    id_z = basis_z.generator(xtype_z).getL2Norm(j_z,k_z);
    id_z *= id_z;
    dd_z = basis_z.generator(xtype_z).getH1SemiNorm(j_z,k_z);
    dd_z *= dd_z;

    return (T)1./std::sqrt(a_x*dd_x*id_y*id_z + a_y*id_x*dd_y*id_z + a_z*id_x*id_y*dd_z
                            + c*id_x*id_y*id_z);
}

template <typename T, typename Basis3D>
T
OptimizedH1Preconditioner3D<T,Basis3D>::operator()(const Index3D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k,
                            index.index3.xtype,index.index3.j,index.index3.k);
}

template <typename T, typename Basis3D>
T
OptimizedH1Preconditioner3D<T,Basis3D>::operator[](const Index3D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k,
                            index.index3.xtype,index.index3.j,index.index3.k);
}

}   // namespace lawa
