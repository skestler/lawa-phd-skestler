/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <cassert>
#include <cmath>
#include <lawa/math/math.h>

namespace lawa {

//--- Gauss-Legendre Quadrature-------------------------------------------------

template <typename Integral>
Quadrature<Gauss,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _order(-1)
{
    _legendre(_precalculatedOrder);
    setOrder(4);
}

template <typename Integral>
const typename Integral::T
Quadrature<Gauss,Integral>::operator()(T a, T b) const
{
    T ret = 0.0;
    for (int i=1; i<=_order; ++i) {
        ret += _weights(_order,i) * integral.integrand(0.5*(b-a)*_knots(_order,i)+0.5*(b+a));
    }
    ret *= 0.5*(b-a);

    return ret;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::setOrder(int order)
{
    assert(order>0);

    if (order>=_precalculatedOrder) {
        _legendre(order);
        _precalculatedOrder = order;
    }
    _order = order;
}

template <typename Integral>
int
Quadrature<Gauss,Integral>::order() const
{
    return _order;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::_legendre(int order)
{
    T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(order, order);
    _weights.engine().resize(order, order);

    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=order; ++k) {
        int     m = (k+1)/2;
        T xm = 0.5 * (x2+x1),
          xl = 0.5 * (x2-x1);
        for (int i=1; i<=m; ++i) {
            T z = cos(M_PI*(i-0.25)/(k+0.5)),
              z1, pp;
            do {
                T p1 = 1.0,
                  p2 = 2.0;
                for (int j=1; j<=k; ++j) {
                    T p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            _knots(k,i)     = xm - xl*z;
            _knots(k,k+1-i) = xm + xl*z;
            _weights(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
            _weights(k,k+1-i) = _weights(k,i);
        }
    }
}

template <typename Integral>
flens::GeMatrix<flens::FullStorage<typename Integral::T,cxxblas::ColMajor> >
Quadrature<Gauss,Integral>::_weights;

template <typename Integral>
flens::GeMatrix<flens::FullStorage<typename Integral::T,cxxblas::ColMajor> >
Quadrature<Gauss,Integral>::_knots;

template <typename Integral>
int
Quadrature<Gauss,Integral>::_precalculatedOrder = 10;

//---  Trapezoidal rule -------------------------------------------------------
template <typename Integral>
Quadrature<Trapezoidal,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _n(-1)
{
    setN(100);
}

template <typename Integral>
const typename Integral::T
Quadrature<Trapezoidal,Integral>::operator()(T a, T b) const
{
    T h = (b-a) / _n;
    T ret = .5 * h * integral.integrand(a);
    a += h;
    for (int i=1; i<_n; ++i, a+=h) {
        ret += h * integral.integrand(a);
    }
    ret += .5 * h * integral.integrand(b);

    return ret;
}

template <typename Integral>
int
Quadrature<Trapezoidal,Integral>::n() const
{
    return _n;
}

template <typename Integral>
void
Quadrature<Trapezoidal,Integral>::setN(const int n)
{
    assert(n>0);

    _n = n;
}

//---  ExpWeighted Rule -------------------------------------------------------
template <typename Integral>
Quadrature<ExpWeighted,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _eta(-0.5*std::log(integral.function(1.))),
      _max_polynomialorder(2)
{
    if (integral._f2) {
        _max_polynomialorder = integral.first.d-1 + integral.second.d-1;
    }
    else {
        _max_polynomialorder = integral.first.d;
    }
    std::cout << "_eta = " << _eta << std::endl;
}

template <typename Integral>
const typename Integral::T
Quadrature<ExpWeighted,Integral>::operator()(T a, T b) const
{
    if (a>=b) return 0.;
    if (_max_polynomialorder<=2) {
        long double p1_a=0., p1_b=0., p2_a=0., p2_b=0., d_p1=0., d_p2=0.;
        long double x = a+0.1*(b-a);
        if (integral.deriv1==0) {
            p1_a = integral.first.generator(integral.e1)
                     ((a-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,0)
                     /(integral.SqrtRightmLeft);
            p1_b = integral.first.generator(integral.e1)
                     ((b-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,0)
                     /(integral.SqrtRightmLeft);
            d_p1 = integral.first.generator(integral.e1)
                    ((x-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,1)
                    / (integral.SqrtRightmLeft * integral.RightmLeft);
        }
        else {
            p1_a= integral.first.generator(integral.e1)
                                ((x-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,1)
                                / (integral.SqrtRightmLeft * integral.RightmLeft);
            p1_b = p1_a;
        }
        if (integral.deriv2==0) {
            p2_a = integral.second.generator(integral.e2)
                     ((a-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,0)
                     /(integral.SqrtRightmLeft);
            p2_b = integral.second.generator(integral.e2)
                     ((b-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,0)
                     /(integral.SqrtRightmLeft);
            d_p2 = integral.second.generator(integral.e2)
                    ((x-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,1)
                    / (integral.SqrtRightmLeft * integral.RightmLeft);
        }
        else {
            p2_a = integral.second.generator(integral.e2)
                    ((x-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,1)
                    / (integral.SqrtRightmLeft * integral.RightmLeft);
            p2_b = p2_a;
        }




        if (b<=0) {
            long double Exp_p_2eta_a = exp((long double)2*_eta*a);
            long double Exp_p_2eta_b = exp((long double)2*_eta*b);

            long double help1 = (p1_b*p2_b*Exp_p_2eta_b-p1_a*p2_a*Exp_p_2eta_a)/((long double)2.*_eta);
            long double help2 = d_p1*(p2_b*Exp_p_2eta_b-p2_a*Exp_p_2eta_a)/((long double)4.*_eta*_eta);
            long double help3 = d_p2*(p1_b*Exp_p_2eta_b-p1_a*Exp_p_2eta_a)/((long double)4.*_eta*_eta);
            long double help4 = 2*d_p1*d_p2*(Exp_p_2eta_b-Exp_p_2eta_a)/((long double)8.*_eta*_eta*_eta);

            return help1-help2-help3+help4;

        }
        else {  // a >= 0
            long double Exp_m_2eta_a = exp((long double)-2*_eta*a);
            long double Exp_m_2eta_b = exp((long double)-2*_eta*b);

            long double help1 = (p1_b*p2_b*Exp_m_2eta_b-p1_a*p2_a*Exp_m_2eta_a)/((long double)-2.*_eta);
            long double help2 = d_p1*(p2_b*Exp_m_2eta_b-p2_a*Exp_m_2eta_a)/((long double)4.*_eta*_eta);
            long double help3 = d_p2*(p1_b*Exp_m_2eta_b-p1_a*Exp_m_2eta_a)/((long double)4.*_eta*_eta);
            long double help4 = 2*d_p1*d_p2*(Exp_m_2eta_b-Exp_m_2eta_a)/((long double)-8.*_eta*_eta*_eta);
            //std::cout << "(" << a << "," << b << "): " << help1 << " " << help2 << " " << help3 << " " << help4 << std::endl;
            return help1-help2-help3+help4;
        }
    }
    else {
        assert(0);
        return 0.;
    }
}

} // namespace lawa

