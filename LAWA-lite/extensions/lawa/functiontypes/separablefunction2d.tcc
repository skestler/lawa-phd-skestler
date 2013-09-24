namespace lawa {

template<typename T>
SeparableFunction2D<T>::SeparableFunction2D(Function<T> _F_x, Function<T> _F_y)
    : F_x(_F_x), F_y(_F_y)
{
}

template<typename T>
SeparableFunction2D<T>::SeparableFunction2D(T (*_f_x)(T),
        const DenseVector<Array<T> > &_singularPts_x,
        T (*_f_y)(T), const DenseVector<Array<T> > &_singularPts_y)
    : F_x(Function<T>(_f_x, _singularPts_x)), F_y(Function<T>(_f_y, _singularPts_y))
{
}

template<typename T>
T
SeparableFunction2D<T>::operator()(T x, T y) const
{
    return F_x(x) * F_y(y);
}

} //  namespace lawa

