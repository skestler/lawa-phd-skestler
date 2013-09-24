/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

namespace lawa {

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
RHS<T,Index,RHSINTEGRAL,Preconditioner>::RHS(const RHSINTEGRAL &_rhsintegral, Preconditioner &_P)
    :    rhsintegral(_rhsintegral), P(_P), rhs_data(), rhs_abs_data()
{
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
RHS<T,Index,RHSINTEGRAL,Preconditioner>::RHS(const RHSINTEGRAL &_rhsintegral, Preconditioner &_P,
                                             const Coefficients<Lexicographical,T,Index> &_rhs_data)
    :    rhsintegral(_rhsintegral), P(_P), rhs_data(), rhs_abs_data()
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    for (const_coeff_it it=_rhs_data.begin(); it!=_rhs_data.end(); ++it) {
        rhs_data[(*it).first] = (*it).second;
        rhs_abs_data.insert(val_type((*it).second, (*it).first));
    }
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
T
RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(const Index &lambda)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;
    const_coeff_it it_end       = rhs_data.end();
    const_coeff_it it_index     = rhs_data.find(lambda);

    if (it_index != it_end) {
        return (*it_index).second;
    }
    else {
        T ret = P(lambda) * rhsintegral(lambda);
        rhs_data[lambda] = ret;
        rhs_abs_data.insert(val_type(ret, lambda));
        return ret;
    }
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(const IndexSet<Index> &Lambda)
{
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    Coefficients<Lexicographical,T,Index> ret;
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(*lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index> ret;
    ret = THRESH(rhs_abs_data,tol,true);
    return ret;
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
T
RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(T t, const Index &lambda)
{
    T ret;
    if (rhs_data.count(lambda) == 0) {
        ret = P(lambda) * rhsintegral(t, lambda);
        rhs_data[lambda] = ret;
    }
    else {
        ret = rhs_data[lambda];
    }
    return ret;
}

template <typename T, typename Index, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(T t, const IndexSet<Index> &Lambda)
{
    typedef typename IndexSet<Index>::iterator const_set_it;
    Coefficients<Lexicographical,T,Index> ret(Lambda.d,Lambda.d_);
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = RHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(t, *lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}



}    //namespace lawa

