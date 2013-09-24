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

template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
APPLY_Helmholtz2D(MA &A, const Coefficients<Lexicographical,T,Index2D> &v, int k)
{
    typedef typename Coefficients<AbsoluteValue,T,Index2D >::const_iterator abs_const_it;
    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator coeff_const_it;
    typedef typename IndexSet<Index2D>::iterator set_const_it;
    typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

    int jmin_x=100,jmin_y=100;
    for (coeff_const_it it=v.begin(); it!=v.end(); ++it) {
        if ((*it).first.index1.xtype == XBSpline) {
            jmin_x = (*it).first.index1.j;
        }
        if ((*it).first.index2.xtype == XBSpline) {
            jmin_y = (*it).first.index2.j;
        }
    }
    //std::cout << "jmin_x = " << jmin_x << ", jmin_y = " << jmin_y << std::endl;
    Coefficients<Lexicographical,T,Index2D > ret(v.d, v.d_);
    if (v.size() > 0) {
        Coefficients<AbsoluteValue,T,Index2D > temp(v.d,v.d_);
        temp = v;

        int s = 0, count = 0;

        for (abs_const_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index2D> Lambda_v(v.d,v.d_);
            Lambda_v=A.c.SparsityPattern((*it).second, jmin_x, jmin_y, k-s, 1, 0);
            for (set_const_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                T prec = A.prec(*mu) * A.prec((*it).second);
                ret[*mu] += prec *
                            A.data_dd_x((*mu).index1, (*it).second.index1) *
                            A.data_id_y((*mu).index2, (*it).second.index2) *
                            (*it).first;
                //ret[*mu] += A(*mu, (*it).second, 1, 0) * (*it).first;
            }
            Lambda_v=A.c.SparsityPattern((*it).second, jmin_x, jmin_y, k-s, 0, 1);
            for (set_const_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                T prec = A.prec(*mu) * A.prec((*it).second);
                ret[*mu] += prec *
                            A.data_id_x((*mu).index1, (*it).second.index1) *
                            A.data_dd_y((*mu).index2, (*it).second.index2) *
                            (*it).first;
                //ret[*mu] += A(*mu, (*it).second, 0, 1) * (*it).first;
            }
            Lambda_v=A.c.SparsityPattern((*it).second, jmin_x, jmin_y, k-s, 0, 0);
            for (set_const_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                T prec = A.prec(*mu) * A.prec((*it).second);
                ret[*mu] += prec *
                            A.data_id_x((*mu).index1, (*it).second.index1) *
                            A.data_id_y((*mu).index2, (*it).second.index2) *
                            (*it).first;
                //ret[*mu] += A(*mu, (*it).second, 0, 0) * (*it).first;
            }
            /*
            Lambda_v=A.c.SparsityPattern((*it).second, jmin_x, jmin_y, k-s, 1, 1);
            for (set_const_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*it).second, 1, 0) * (*it).first
                           +A(*mu, (*it).second, 0, 1) * (*it).first
                           +A(*mu, (*it).second, 0, 0) * (*it).first;
            }
            */
            ++count;
            //std::cout << (*it).second << ", (" << count << ", " << v.size() << ") : " << Lambda_v.size() << std::endl;
            s = int(log(T(count))/log(T(2))) + 1;
        }
    }
    return ret;
}


}    //namespace lawa

