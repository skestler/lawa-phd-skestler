namespace lawa {

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
         typename S2_x, typename S2_y>
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::UniformTensorMatrix2D
                                                       (const UniformBasis2D &basis,
                                                        const S1_x &s1_x, const S1_y &s1_y,
                                                        const S2_x &s2_x, const S2_y &s2_y,
                                                        const int Jx, const int Jy)
    : _basis(basis),
      _s1_x(s1_x), _s1_y(s1_y), _s2_x(s2_x), _s2_y(s2_y),
      _Jx(Jx), _Jy(Jy),
      assembler_x(basis.first), assembler_y(basis.second)
{
    this->assembleMatrices();
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
       typename S2_x, typename S2_y>
int
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::numRows() const
{
    return _basis.second.mra.cardI(_Jy)*_basis.first.mra.cardI(_Jx);
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
       typename S2_x, typename S2_y>
int
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::numCols() const
{
    return _basis.second.mra.cardI(_Jy)*_basis.first.mra.cardI(_Jx);
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
       typename S2_x, typename S2_y>
IndexSet<Index2D>
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::getIndexSet() const
{
    IndexSet<Index2D> Lambda;

    for (int k_x=_basis.first.mra.rangeI(_basis.first.j0).firstIndex(); k_x<=_basis.first.mra.rangeI(_basis.first.j0).lastIndex(); ++k_x) {
        for (int k_y=_basis.second.mra.rangeI(_basis.second.j0).firstIndex(); k_y<=_basis.second.mra.rangeI(_basis.second.j0).lastIndex(); ++k_y) {
            Lambda.insert(Index2D(Index1D(_basis.first.j0,k_x,XBSpline),Index1D(_basis.second.j0,k_y,XBSpline)));
        }
        for (int j_y=_basis.second.j0; j_y<=_Jy-1; ++j_y) {
            for (int k_y=_basis.second.rangeJ(j_y).firstIndex(); k_y<=_basis.second.rangeJ(j_y).lastIndex(); ++k_y) {
                Lambda.insert(Index2D(Index1D(_basis.first.j0,k_x,XBSpline),Index1D(j_y,k_y,XWavelet)));
            }
        }
    }
    for (int j_x=_basis.first.j0; j_x<=_Jx-1; ++j_x) {
        for (int k_x=_basis.first.rangeJ(j_x).firstIndex(); k_x<=_basis.first.rangeJ(j_x).lastIndex(); ++k_x) {
            for (int k_y=_basis.second.mra.rangeI(_basis.second.j0).firstIndex(); k_y<=_basis.second.mra.rangeI(_basis.second.j0).lastIndex(); ++k_y) {
                Lambda.insert(Index2D(Index1D(j_x,k_x,XWavelet),Index1D(_basis.second.j0,k_y,XBSpline)));
            }
            for (int j_y=_basis.second.j0; j_y<=_Jy-1; ++j_y) {
                for (int k_y=_basis.second.rangeJ(j_y).firstIndex(); k_y<=_basis.second.rangeJ(j_y).lastIndex(); ++k_y) {
                    Lambda.insert(Index2D(Index1D(j_x,k_x,XWavelet),Index1D(j_y,k_y,XWavelet)));
                }
            }
        }
    }
    return Lambda;
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
        typename S2_x, typename S2_y>
flens::DenseVector<flens::Array<T> >
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::operator*(const DenseVectorT &v) const
{
    //Timer time;
    //time.start();
    if (v.length()!=_basis.second.mra.cardI(_Jy)*_basis.first.mra.cardI(_Jx)) {
        std::cerr << "UniformTensorMatrix2D: Dimension mismatch!" << std::endl;
        exit(1);
    }

    DenseMatrixT V(_basis.second.mra.cardI(_Jy), _basis.first.mra.cardI(_Jx));
    for (int j=1; j<=_basis.first.mra.cardI(_Jx); ++j) {
        int firstIndex=(j-1)*_basis.second.mra.cardI(_Jy)+1;
        int lastIndex =    j*_basis.second.mra.cardI(_Jy);
        V(_,j) = v(_(firstIndex,lastIndex));
    }

    DenseMatrixT Ret = this->operator*(V);

    DenseVectorT ret(_basis.second.mra.cardI(_Jy) * _basis.first.mra.cardI(_Jx));
    for (int j=1; j<=_basis.first.mra.cardI(_Jx); ++j) {
        int firstIndex=(j-1)*_basis.second.mra.cardI(_Jy)+1;
        int lastIndex =    j*_basis.second.mra.cardI(_Jy);
        ret(_(firstIndex,lastIndex)) = Ret(_,j);
    }
    //time.stop();
    //std::cerr << "DEBUG: Time for tensor matrix vector multiplication " << time.elapsed() << std::endl;
    return ret;
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
        typename S2_x, typename S2_y>
void
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::assembleMatrices()
{
    _M_s1_x = assembler_x.assembleStiffnessMatrix(_s1_x, _Jx);
    _M_s1_y   = assembler_y.assembleStiffnessMatrix(_s1_y, _Jy);
    _M_s2_x = assembler_x.assembleStiffnessMatrix(_s2_x, _Jx);
    _M_s2_y   = assembler_y.assembleStiffnessMatrix(_s2_y, _Jy);
}

template<typename T, typename UniformBasis2D, typename S1_x, typename S1_y,
        typename S2_x, typename S2_y>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
UniformTensorMatrix2D<T,UniformBasis2D,S1_x,S1_y,S2_x,S2_y>::operator*(const DenseMatrixT &V) const
{
    DenseMatrixT ret(_basis.second.mra.cardI(_Jy), _basis.first.mra.cardI(_Jx));
    if ((V.numRows()!=_basis.second.mra.cardI(_Jy)) || (V.numCols()!=_basis.second.mra.cardI(_Jx)))
    {
        std::cerr << "UniformTensorMatrix2D: Dimension mismatch!" << std::endl;
        exit(1);
    }
    for (int j=1; j<=V.numCols(); ++j) {
        DenseVectorT tmp = _M_s1_y*V(_,j);
        ret(_,j) = tmp;
    }
    for (int i=1; i<=ret.numRows(); ++i) {
        DenseVectorT v=ret(i,_);
        DenseVectorT tmp = _M_s1_x*v;
        ret(i,_) = tmp;
    }


    DenseMatrixT intermed(_basis.second.mra.cardI(_Jy), _basis.first.mra.cardI(_Jx));
    for (int j=1; j<=V.numCols(); ++j) {
        DenseVectorT tmp = _M_s2_y*V(_,j);
        intermed(_,j) = tmp;
    }

    for (int i=1; i<=intermed.numRows(); ++i) {
        DenseVectorT v=intermed(i,_);
        DenseVectorT tmp = _M_s2_x*v;
        ret(i,_) += tmp;
    }
    return ret;
}

}   // namespace lawa
