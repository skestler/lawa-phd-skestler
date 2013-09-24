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

namespace lawa {
    
template <typename MA, typename VX, typename VB>
int
gmres(const MA &A, VX &x, const VB &b, typename _gmres<MA>::T tol,
      long maxIterations)
{
    typedef typename _gmres<MA>::T              T;
    typedef DenseVector<Array<T> >              DeVector;
    typedef GeMatrix<FullStorage<T, cxxblas::ColMajor> > DeMatrix;

    int N = b.length();
    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    DeVector r = b - A*x;
    T beta = sqrt(r*r);
    if (beta==0.0) {
        return 0;
    }

    DeMatrix  V(N, maxIterations+1);
    //V(_,1) = r / beta;
    V(_,1) = 1./beta * r;
    
    DeMatrix  H(maxIterations+1, maxIterations  );

    DeVector w_j(N);
    DeVector g(maxIterations+1);
    DeVector c(maxIterations+1), s(maxIterations+1);

    T nu, rho = tol + 1;
    T Htemp, h_ij;
    
    g(1) = beta;

    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }

        w_j =  A * V(_,j);
        T normInitialWj = sqrt(w_j * w_j);

        for (int i=1; i<=j; ++i) {
            H(i,j) = w_j * V(_,i);
            w_j = w_j - H(i,j) * V(_, i);
        }
          H(j+1,j) = sqrt(w_j*w_j);

        if (H(j+1,j) / normInitialWj < 1.0) {
            for (int i=1; i<=j; ++i) {
                Htemp = w_j * V(_, i) ;
                w_j = w_j - Htemp * V(_, i);
            }
            H(j+1, j) = sqrt(w_j * w_j);
          }


        for (int i=1; i<=j-1; ++i) {
            h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
            H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
            H(i,j) =    h_ij;
        }

        nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
        if (nu!=0.0) {
            //V(_,j+1) = w_j / H(j+1,j);
            V(_,j+1) = 1. / H(j+1,j) * w_j;
        }

        if (nu!=0.0) {
            s(j) = H(j+1,j)/nu;
            c(j) = H(j,j)/nu;
            H(j,j) = nu;
            H(j+1,j) = 0.0;
            g(j+1) = -s(j)*g(j);
            g(j) = c(j)*g(j);
        }
        rho = fabs(g(j+1));
    }

    if (j>=1) {
        DeVector y(j);
        for (int i=j; i>=1; --i) {
            y(i) = g(i) / H(i,i);
            for (int l=j; l>i; --l) {
                y(i) -= H(i,l) * y(l) / H(i,i);
            }
          }

        //todo: For some unknown reasons, this causes trouble in waveletgalerkinoptionpricer when
        //      long double and a sufficiently large number of time steps is used...

        //x += V(_,_(1,j)) * y;

        DeMatrix tmpV(N,j);
        tmpV = V(_,_(1,j));

        DeVector temp(N);
        for (int i=1; i<=N; ++i) {
            for (int k=1; k<=j; ++k) {
                temp(i) += tmpV(i,k) * y(k);
            }
        }
        x += temp;

    }
    return j;
}

template <typename MA, typename VX, typename VB>
int
gmresm(const MA &A, VX &x, const VB &b, typename _gmres<MA>::T tol,
       int restart, long maxIterations)
{
    typedef typename _gmres<MA>::T T;
    typedef DenseVector<Array<T> > DeVector;

    int N=b.length();
    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
    DeVector r=b-A*x;
    while (sqrt(r*r)>tol && k<=maxIterations) {
        k+=gmres(A,x,b,tol,restart);
        r=A*x-b;
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << sqrt(r*r) << std::endl;
        #endif
    }
    return k;
}

template <typename Prec, typename MA, typename VX, typename VB>
int pgmres (const Prec &P, const MA &A, VX &x, const VB &b,
            typename _gmres<MA>::T tol, long maxIterations)
{
    typedef typename _gmres<MA>::T              T;
    typedef DenseVector<Array<T> >              DeVector;
    typedef GeMatrix<FullStorage<T, cxxblas::ColMajor> > DeMatrix;

    int N = b.length();
    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    DeVector temp = b - A*x, r = transpose(P)*temp;
    T beta = sqrt(r * r);
    if (beta==0.0) {
        return 0;
    }

    DeMatrix  V(N,maxIterations + 1);
    //V(_,1) = r / beta;
    V(_,1) = 1. / beta * r;
    DeMatrix  H(maxIterations + 1,maxIterations    );

    DeVector w_j(N);
    DeVector g(maxIterations + 1);
    DeVector c(maxIterations + 1), s(maxIterations + 1);

    T nu, rho = tol + 1;
    T Htemp, h_ij;

    g(1) = beta;

    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }
          
          w_j = P*V(_,j);
        temp = A*w_j;
        w_j = transpose(P)*temp;
        
        T normInitialWj = sqrt(w_j * w_j);

        for (int i=1; i<=j; ++i) {
            H(i, j) = w_j * V(_,i);
            w_j = w_j - H(i, j) * V(_, i);
        }
          H(j+1, j) = sqrt(w_j * w_j);

        if ( H(j+1, j) / normInitialWj < 1.0) {
            for (int i=1; i<=j; ++i) {
                Htemp = w_j * V(_, i) ;
                w_j = w_j - Htemp * V(_, i);
            }
            H(j+1, j) = sqrt(w_j*w_j);
          }

        for (int i=1; i<=j-1; ++i) {
            h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
            H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
            H(i,j) =    h_ij;
        }

        nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
        if (nu!=0.0) {
            //V(_,j+1) = w_j / H(j+1,j);
            V(_,j+1) = 1. / H(j+1,j) * w_j;
        }

        if (nu!=0.0) {
            s(j) = H(j+1,j)/nu;
            c(j) = H(j,j)/nu;
            H(j,j) = nu;
            H(j+1,j) = 0.0;
            g(j+1) = -s(j)*g(j);
            g(j) = c(j)*g(j);
        }
        rho = fabs(g(j+1));
    }

    if (j >= 1 ) {
        DeVector y(j);
        for (int i=j; i>=1; --i) {
            y(i) = g(i) / H(i,i);
            for (int l=j; l>i; --l) {
                y(i) -= H(i,l) * y(l) / H(i,i);
            }
          }

        //todo: For some unknown reasons, this causes trouble in waveletgalerkinoptionpricer when
        //      long double and a sufficiently large number of time steps is used...
        //temp = V(_,_(1,j)) * y;

        DeMatrix tmpV(N,j);
        tmpV = V(_,_(1,j));

        for (int i=1; i<=N; ++i) {
            temp(i) = 0.;
            for (int k=1; k<=j; ++k) {
                temp(i) += tmpV(i,k) * y(k);
            }
        }
        x += P * temp;
    }
    return j;
}

template <typename Prec, typename MA, typename VX, typename VB>
int
pgmresm(const Prec &Pr, const MA &A, VX &x, const VB &b,
        typename _gmres<MA>::T tol, int restart, long maxIterations)
{
    typedef typename _gmres<MA>::T T;
    typedef DenseVector<Array<T> > DeVector;

    int N=b.length();
    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
    DeVector r=b-A*x;
    while ((sqrt(r*r)>tol) && (k<=maxIterations)) {
        k+=pgmres(Pr,A,x,b,tol,restart);
        r=A*x-b;
    }
    return k;
}

} // namespace lawa

