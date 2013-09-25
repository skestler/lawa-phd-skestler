#ifndef _ND_GENERATOR2D_HPP
#define _ND_GENERATOR2D_HPP 1

#include "cholesky.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>


using namespace boost::numeric::ublas;


template <typename T>
class ND_Generator2D
{
	
	typedef boost::variate_generator<boost::mt19937, boost::normal_distribution<> > ND_Generator1D;
	
public:
    ND_Generator2D(matrix<T> _Q, vector<T> _drift)
    : Q(_Q), drift(_drift), nd_generator1d(boost::mt19937(), boost::normal_distribution<>())
    {
        L.resize(Q.size1(), Q.size2());
        size_t n = cholesky_decompose(Q, L);
        //std::cerr << "Q = " << Q << ", L = " << L << std::endl;
    }

   ND_Generator2D(T q11, T q12, T q21, T q22, T drift1, T drift2)
    : nd_generator1d(boost::mt19937(), boost::normal_distribution<>())
    {
        Q.resize(2,2);
        Q(0,0) = q11; Q(0,1) = q12; Q(1,0) = q21; Q(1,1) = q22;
        L.resize(2,2);

        drift.resize(2);
        drift(0) = drift1;
        drift(1) = drift2;

        if (q12!=q21) {
            std::cerr << "ND_Generator2D: Q is not symmetric!" << std::endl;
            exit(1);
        }
        if (fabs(q21)<1e-12) {
            std::cerr << "ND_Generator2D: Q is already a diagonal matrix!" << std::endl;
            L(0,0) = sqrt(q11);
            L(1,1) = sqrt(q22);
            L(0,1) = 0.;
            L(1,0) = 0.;
        }
        else {
            size_t n = cholesky_decompose(Q, L);
            L(0,1) = 0.;
        }
        std::cerr.precision(16);
        std::cerr << "Q = " << Q << ", L = " << L << std::endl;
    }
	
	T
	operator()(T &y1, T &y2)
	{
		T x1 = nd_generator1d();
		T x2 = nd_generator1d();
		
		y1 = drift(0) + L(0,0)*x1 + L(0,1)*x2; 
		y2 = drift(1) + L(1,0)*x1 + L(1,1)*x2; 
	}
	
	T
	operator()(T maturity, T &y1, T &y2)
	{
		T x1 = nd_generator1d();
		T x2 = nd_generator1d();
		y1 = drift(0)*maturity + std::sqrt(maturity)*(L(0,0)*x1 + L(0,1)*x2); 
		y2 = drift(1)*maturity + std::sqrt(maturity)*(L(1,0)*x1 + L(1,1)*x2); 
	}
	
	matrix<T> Q;
	matrix<T> L;
	vector<T> drift;
	ND_Generator1D nd_generator1d;
};

#endif
