#ifndef BENCH_BENCHMARK_H
#define BENCH_BENCHMARK_H 1

template <typename Operation>
struct Benchmark
{
    Benchmark(long _minProblemSize = 20,  long _maxProblemSize = 10000,
              long _numberOfMeasurements = 100,
              long _minNumberOfComputations = 30,
              double _minTimeToElapse = 2.0);

    void
    run();

    void
    runTests();

    long
    _computeProblemSize(long i);

    long minProblemSize;
    long maxProblemSize;
    long numberOfMeasurements;
    long minNumberOfComputations;
    double minTimeToElapse;
    Timer           timer;
};

#include <bench/benchmark.tcc>

#endif // BENCH_BENCHMARK_H
