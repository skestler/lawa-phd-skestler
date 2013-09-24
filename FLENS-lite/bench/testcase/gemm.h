#ifndef BENCH_TESTCASE_GEMM_H
#define BENCH_TESTCASE_GEMM_H 1

template <typename M>
struct Gemm
{
    typedef M                               Matrix;
    typedef typename Matrix::IndexType      IndexType;
    typedef typename Matrix::ElementType    ElementType;

    Gemm(long _N);

    void
    run(long numberOfComputations);

    long
    numBaseOperations();

    bool
    test();

    //--> the following two method need to be implemeneted for each benchmark
    static const char *
    name();

    void
    compute();
    //<--

    long            N;
    Matrix          A, B, C;
    ElementType     alpha, beta;
};

#include <bench/testcase/gemm.tcc>

#endif // BENCH_TESTCASE_GEMM_H
