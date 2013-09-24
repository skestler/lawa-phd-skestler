namespace lawa{

template<typename T, typename Basis>
BlockAssembler1D<T, Basis>::BlockAssembler1D(const Basis& _basis)
    : basis(_basis)
{
}

template<typename T, typename Basis>
template <typename BilinearForm>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
BlockAssembler1D<T, Basis>::assembleStiffnessMatrixBlock(BilinearForm& a, int i1, int i2, T tol)
{
    assert(i1>=-1);
    assert(i2>=-1);
    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;

    int N1=0, N2=0;
    if (i1<0) { N1 = basis.mra.cardI(j0); }
    else      { N1 = basis.cardJ(j0+i1); }
    if (i2<0) { N2 = basis.mra.cardI(j0); }
    else      { N2 = basis.cardJ(j0+i2); }

    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A(N1,N2);

    if (i1<0) {
        if (i2<0) {    //SF*SF
            for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
                for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
                    T val = a(XBSpline, j0, k1, XBSpline, j0, k2);
                    if(fabs(val) > tol){
                        A(k1-offsetI, k2-offsetI) = val;
                    }
                }
            }
        }

        else {    // SF * W
            for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
                for(int k2 = basis.rangeJ(j0+i2).firstIndex(); k2 <= basis.rangeJ(j0+i2).lastIndex(); ++k2){
                    T val = a(XBSpline, j0, k1, XWavelet, j0+i2, k2);
                    if(fabs(val) > tol){
                        A(k1-offsetI,  k2-offsetJ) = val;
                    }
                }
            }
        }
    }

    else {
        if (i2<0) { // W * SF
            for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
                for(int k1 = basis.rangeJ(j0+i1).firstIndex(); k1 <= basis.rangeJ(j0+i1).lastIndex(); ++k1){
                    T val = a(XWavelet, j0+i1, k1, XBSpline, j0, k2);
                    if(fabs(val) > tol){
                        A(k1-offsetJ, k2-offsetI) = val;
                    }
                }
            }
        }
        else {    // W * W
            for(int k1 = basis.rangeJ(j0+i1).firstIndex(); k1 <= basis.rangeJ(j0+i1).lastIndex(); ++k1){
                for(int k2 = basis.rangeJ(j0+i2).firstIndex(); k2 <= basis.rangeJ(j0+i2).lastIndex(); ++k2){
                    T val = a(XWavelet, j0+i1, k1, XWavelet, j0+i2, k2);
                    if(fabs(val) > tol){
                        A(k1-offsetJ, k2-offsetJ) = val;
                    }
                }
            }
        }
    }

    A.finalize();

    return A;
}

template<typename T, typename Basis>
template <typename BilinearForm>
flens::DenseVector<flens::Array<T> >
BlockAssembler1D<T, Basis>::assembleStiffnessMatrixBlockDiagonal(BilinearForm& a, int i, T tol)
{
    assert(i>=-1);

    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;

    int N=0;
    if (i<0) { N = basis.mra.cardI(j0); }
    else     { N = basis.cardJ(j0+i); }


    flens::DenseVector<flens::Array<T> > diag_A(N);

    if (i<0) {
        for(int k = basis.mra.rangeI(j0).firstIndex(); k <= basis.mra.rangeI(j0).lastIndex(); ++k){
            T val = a(XBSpline, j0, k, XBSpline, j0, k);
            if(fabs(val) > tol){
                diag_A(k-offsetI) = val;
            }
        }
    }

    else {
        for(int k = basis.rangeJ(j0+i).firstIndex(); k <= basis.rangeJ(j0+i).lastIndex(); ++k){
            T val = a(XWavelet, j0+i, k, XWavelet, j0+i, k);
            if(fabs(val) > tol){
                diag_A(k-offsetJ) = val;
            }
        }
    }

    return diag_A;
}

template<typename T, typename Basis>
template <typename RHSIntegral>
flens::DenseVector<flens::Array<T> >
BlockAssembler1D<T, Basis>::assembleRHSBlock(RHSIntegral &rhs, int i)
{
    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;

    int N=0;
    if (i<0) { N = basis.mra.cardI(j0); }
    else     { N = basis.cardJ(j0+i); }

    flens::DenseVector<flens::Array<T> > f(N);
    if (i<0) {
        for(int k = basis.mra.rangeI(j0).firstIndex(); k <= basis.mra.rangeI(j0).lastIndex(); ++k){
            f(k - offsetI) = rhs(XBSpline, j0, k);
        }
    }
    else {
        for(int k = basis.rangeJ(j0+i).firstIndex(); k <= basis.rangeJ(j0+i).lastIndex(); ++k){
            f(k - offsetJ) = rhs(XWavelet, j0+i, k);
        }
    }
    return f;
}



}

