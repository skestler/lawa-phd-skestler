/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_ADAPTIVEPAYOFFQUADRATURE2D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_ADAPTIVEPAYOFFQUADRATURE2D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <applications/finance/options/options.h>
#include <applications/finance/processes/processes.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

namespace lawa {

struct my_params
{
    double a, b, c, d;
};

double
sing_basketput2d_f (double x, void * params)
{
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    return exp (b*x) + d*exp (c*x) - a;
}

double
sing_basketput2d_df (double x, void * params)
{
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    return b*exp (b*x) + c*d*exp (c*x);
}

void
sing_basketput2d_fdf (double x, void * params, double * f, double * df)
{
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    double tmp1 = exp (b*x);
    double tmp2 = d*exp (c*x);

    *f  = tmp1 + tmp2 - a;
    *df = b*tmp1 + c*tmp2;   /* uses existing value */
}


template <OptionTypenD OType, typename PayoffIntegral>
struct AdaptivePayoffQuadrature2D
{

};

template <typename PayoffIntegral>
struct AdaptivePayoffQuadrature2D<BasketPut,PayoffIntegral>
{
    typedef typename PayoffIntegral::T T;

    AdaptivePayoffQuadrature2D(PayoffIntegral &_payoffintegral);

    T
    integrate(T a1, T b1, T a2, T b2);

    T
    integrate_smooth(T a1, T b1, T a2, T b2);

    T
    integrate_singular(T a1, T b1, T a2, T b2);

    bool
    contains_singularity(T a1, T b1, T a2, T b2);

    T
    find_intersectionpoint_y1_given_y2(T y2);

    T
    find_intersectionpoint_y2_given_y1(T y1);

    PayoffIntegral &payoffintegral;

    T thresh;

    const gsl_root_fdfsolver_type * rootSolverType;
    gsl_root_fdfsolver            * rootSolver;
    gsl_function_fdf              FDF;
};

}   //namespace lawa

#include <applications/finance/initialconditions/adaptivepayoffquadrature2d.tcc>

#endif  // APPLICATIONS_FINANCE_INITIALCONDITIONS_ADAPTIVEPAYOFFQUADRATURE2D_H
