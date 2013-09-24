#include <bench/aux/random.h>

namespace flens {

template <typename E>
    struct StorageInfo;

} // namespace flens

//------------------------------------------------------------------------------

template <typename M>
Gemm<M>::Gemm(long _N)
    : N(_N), A(N, N), B(N, N), C(N, N)
{
    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            A(i, j) = randomValue<ElementType>();
            B(i, j) = randomValue<ElementType>();
            C(i, j) = randomValue<ElementType>();
        }
    }
    alpha = 1;
    beta = 0;
}

template <typename M>
void
Gemm<M>::run(long numberOfComputations)
{
    for (long i=0; i<numberOfComputations; ++i) {
        compute();
    }
}

template <typename M>
long
Gemm<M>::numBaseOperations()
{
    return 2*N*N*N;
}

template <typename M>
bool
Gemm<M>::test()
{
#ifdef HAVE_CBLAS
    Matrix          ARef = A, BRef = B, CRef = C;
    ElementType     alphaRef = alpha;
    ElementType     betaRef = beta;
    ElementType     tol = std::numeric_limits<ElementType>::epsilon();

    cxxblas::gemm(flens::StorageInfo<typename Matrix::Engine>::Order,
                  cxxblas::NoTrans, cxxblas::NoTrans,
                  CRef.numRows(), CRef.numCols(), ARef.numCols(),
                  alphaRef,
                  ARef.engine().data(), ARef.engine().leadingDimension(),
                  BRef.engine().data(), BRef.engine().leadingDimension(),
                  betaRef,
                  CRef.engine().data(), CRef.engine().leadingDimension());

    compute();

    for (IndexType i=C.firstRow(); i<=C.lastRow(); ++i) {
        for (IndexType j=C.firstCol(); j<=C.lastRow(); ++j) {
            if (abs(C(i,j)-CRef(i,j))>tol) {
                return false;
            }
        }
    }
    return true;
#else 
    std::cerr << "link with native BLAS for testing" << std::endl;
    return false;
#endif
}
