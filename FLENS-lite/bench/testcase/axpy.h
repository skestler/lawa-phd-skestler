#ifndef BENCH_TESTCASE_AXPY_H
#define BENCH_TESTCASE_AXPY_H 1

template <typename V>
struct Axpy
{
    typedef V                               Vector;
    typedef typename Vector::IndexType      IndexType;
    typedef typename Vector::ElementType    ElementType;

    Axpy(long _N);

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
    Vector          x, y;
    ElementType     alpha;
};

#include <bench/testcase/axpy.tcc>

#endif // BENCH_TESTCASE_AXPY_H
