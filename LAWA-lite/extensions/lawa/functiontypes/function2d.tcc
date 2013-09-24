namespace lawa {

template <typename T>
Function2D<T>::Function2D(T (*_f)(T,T), const DenseVector<Array<T> > &_singularPts_x,
         const DenseVector<Array<T> > &_singularPts_y)
    : f(_f), singularPts_x(_singularPts_x), singularPts_y(_singularPts_y)
{
}

template <typename T>
T
Function2D<T>::operator()(T x, T y) const
{
    return f(x,y);
}

}   //namespace lawa

