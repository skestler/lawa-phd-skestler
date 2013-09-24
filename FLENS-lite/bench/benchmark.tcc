#include <cmath>
#include <iostream>

template <typename Operation>
Benchmark<Operation>::Benchmark(long _minProblemSize,  long _maxProblemSize,
                                long _numberOfMeasurements,
                                long _minNumberOfComputations,
                                double _minTimeToElapse)
    : minProblemSize(_minProblemSize),
      maxProblemSize(_maxProblemSize),
      numberOfMeasurements(_numberOfMeasurements),
      minNumberOfComputations(_minNumberOfComputations),
      minTimeToElapse(_minTimeToElapse)
{
}

template <typename Operation>
void
Benchmark<Operation>::run()
{
    std::cout << Operation::name() << std::endl;

    for (long i=0; i<numberOfMeasurements; ++i) {
        long numberOfComputations = minNumberOfComputations;
        long N = _computeProblemSize(i);

        while (true) {
            Operation  op(N);

            timer.start();
            op.run(numberOfComputations);
            timer.stop();

            double elapsed = timer.elapsed();

            if (elapsed<minTimeToElapse) {
                numberOfComputations *= 2;
                continue;
            }

            double mflops = numberOfComputations*op.numBaseOperations();
            mflops /= elapsed*1000000;

            std::cout << i << " "
                      << N << " "
                      << numberOfComputations << " "
                      << elapsed << " "
                      << mflops
                      << std::endl;
            break;
        }
    }
}

template <typename Operation>
void
Benchmark<Operation>::runTests()
{
    std::cout << Operation::name() << std::endl;

    for (long i=0; i<numberOfMeasurements; ++i) {
        long N = _computeProblemSize(i);
        Operation  op(N);

        std::cout << i << ") " << "N = " << N << " ... ";

        if (! op.test()) {
            std::cout << "failed." << std::endl;
            break;
        }
        std::cout << "passed." << std::endl;
    }
}

template <typename Operation>
long
Benchmark<Operation>::_computeProblemSize(long i)
{
    const double ratioProblemSize = double(maxProblemSize)/minProblemSize;
    const double exponent = double(i)/(numberOfMeasurements-1);

    return minProblemSize*pow(ratioProblemSize, exponent);
}
