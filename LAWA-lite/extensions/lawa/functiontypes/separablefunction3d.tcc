namespace lawa {

template<typename T>
SeparableFunction3D<T>::SeparableFunction3D(Function<T> _F_x, Function<T> _F_y, Function<T> _F_z)
    : F_x(_F_x), F_y(_F_y), F_z(_F_z)
{
}

template<typename T>
SeparableFunction3D<T>::SeparableFunction3D(T (*_f_x)(T), const DenseVector<Array<T> > &_singularPts_x,
                                            T (*_f_y)(T), const DenseVector<Array<T> > &_singularPts_y,
                                            T (*_f_z)(T), const DenseVector<Array<T> > &_singularPts_z)
    : F_x(Function<T>(_f_x, _singularPts_x)), F_y(Function<T>(_f_y, _singularPts_y)), F_z(Function<T>(_f_z, _singularPts_z))
{
}

template<typename T>
T
SeparableFunction3D<T>::operator()(T x, T y, T z) const
{
    return F_x(x) * F_y(y) * F_z(z);
}


} //  namespace lawa

