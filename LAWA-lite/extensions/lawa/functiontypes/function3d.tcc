namespace lawa {

template <typename T>
Function3D<T>::Function3D(T (*_f)(T,T,T), const DenseVector<Array<T> > &_singularPts_x,
         const DenseVector<Array<T> > &_singularPts_y, const DenseVector<Array<T> > &_singularPts_z)
    : f(_f), singularPts_x(_singularPts_x), singularPts_y(_singularPts_y), singularPts_z(_singularPts_z)
{
}

template <typename T>
T
Function3D<T>::operator()(T x, T y, T z) const
{
    return f(x,y,z);
}

}   //namespace lawa

