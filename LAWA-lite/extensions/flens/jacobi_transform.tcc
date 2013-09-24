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

#include <cmath>
#include <iostream>
#include <vector>

namespace flens {

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
    a[k][l]=h+s*(g-h*tau);

template <typename T>
void
jacobi(T **a, int n, T d[], T **v, int *nrot)
{
    int j, iq, ip,i;
    T tresh, theta, tau, t, sm, s, h, g, c;
    
    std::vector<T> b(n+1);   // only 1..n needed!!!
    std::vector<T> z(n+1);   // only 1..n needed!!!
    for (ip=1; ip<=n; ++ip) {
        for (iq=1; iq<=n; ++iq) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.;
    }
    for (ip=1; ip<=n; ++ip) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    
    for (i=1; i<=50; ++i) {
        sm = 0.0;
        for (ip=1; ip<=n-1; ++ip) {
            for (iq=ip+1; iq<=n; ++iq) {
                sm += fabs(a[ip][iq]);
            }
        }
        if (sm==0.0) {
            return;
        }
        if (i<4) {
            tresh = 0.2*sm/(n*n);
        } else {
            tresh = 0.0;
        }
        for (ip=1; ip<=n-1; ++ip) {
            for (iq=ip+1; iq<=n; ++iq) {
                g = 100.0*fabs(a[ip][iq]);
                if (i>4 && (T)(fabs(d[ip])+g)==(T)fabs(d[ip])
                        && (T)(fabs(d[iq])+g)==(T)fabs(d[iq])) {
                    a[ip][iq] = 0.0;
                } else if (fabs(a[ip][iq])>tresh) {
                    h = d[iq]-d[ip];
                    if ((T)(fabs(h)+g)==(T)fabs(h)) {
                        t = (a[ip][iq])/h;
                    } else {
                        theta = 0.5*h/(a[ip][iq]);
                        t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta<0.0) {
                            t = -t;
                        }
                    }
                    c = 1./sqrt(1+t*t);
                    s = t*c;
                    tau = s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j=1; j<=ip-1; ++j) {
                        ROTATE(a,j,ip,j,iq);
                    }
                    for (j=ip+1; j<=iq-1; ++j) {
                        ROTATE(a,ip,j,j,iq);
                    }
                    for (j=iq+1; j<=n; ++j) {
                        ROTATE(a,ip,j,iq,j);
                    }
                    for (j=1; j<=n; ++j) {
                        ROTATE(v,j,ip,j,iq);
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip=1; ip<=n; ++ip) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    std::cerr << "Too many iterations in routine jacobi" << std::endl;
}

} // namespace flens

