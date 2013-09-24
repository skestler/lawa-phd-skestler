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

#include <fstream>
#include <sstream>

namespace lawa {

template <typename T>
void
spy(const SparseGeMatrix<CRS<T,CRS_General> > &A, const char* filename,
    bool absoluteValues, T eps)
{
    using namespace flens;
    using namespace std;

    stringstream spy_filename;
    spy_filename << filename << ".gps";
    ofstream gps(spy_filename.str().c_str());
    stringstream content;

    //const CRS<T> crs=A.engine();
    int n=A.engine().numCols(), m=A.engine().numRows();
    // get first entry
    T min, max;
    if (absoluteValues==true) {
        min = fabs(A.engine().values(1)); max = fabs(A.engine().values(1));
    } else {
        min = A.engine().values(1); max = A.engine().values(1);
    }

    int nnz = 0;
    T value = 0;
    for (int i=1; i<=m; ++i) {
        for (int j=A.engine().rows(i); j<A.engine().rows(i+1); ++j) {
            if (absoluteValues==true) {
                value = fabs(A.engine().values(j));
            } else {
                value = A.engine().values(j);
            }

            if (fabs(value) >= eps) {
                ++nnz;
                if ( value > max ) {
                    max = value;
                } else if ( value < min ) {
                    min = value;
                }
            }
            // i=column => y-value, j=row => x-value
            content << A.engine().columns(j) << " " << i << " " << value << endl;
        }
    }
    
    gps << "reset; set term png; set output '" << filename << ".png'" << endl;
    if (n==m) {
        gps << "set size square" << endl;
        gps << "set title '$n = " << n << ", n^2 = " << n*n;
    } else {
        gps << "set size ratio " << T(m)/T(n) << endl;
        gps << "set title '$m = " << m << ", n = " << n << ", m*n = " << n*m;
    }

    gps << ", nnz = " << nnz << "$'" << endl
        << "set xrange [1:" << n << "]" << endl
        << "set yrange [1:" << m << "] reverse" << endl
        << "unset xtics; unset ytics; unset ztics" << endl
        << "set pm3d map" << endl
        << "set cbrange[" << min << ":" << max << "]" << endl;

    gps << "set palette model RGB" << endl;
    if (absoluteValues==true) {
        T cbtics = exp(1.0/6.0 * log(max/min));
        gps << "set logscale cb" << endl;
        gps << "set cbtics " << min << ", " << cbtics << ", " << max << endl;
    } else {
        T cbtics = (max-min)/6.0;
        gps << "set cbtics " << min << ", " << cbtics << ", " << max << endl;
    }
    gps << "splot '-' with points pointtype 7 pointsize 0.5 palette notitle"
        << endl << content.str().c_str() << endl;
    gps << "set output" << endl
        << "set term pop" << endl;

    gps.close();
}

template <typename I>
void
spy(const Matrix<I> &A, const char* filename, bool absoluteValues,
    typename I::ElementType eps)
{
    using namespace flens;
    using namespace std;
    typedef typename I::ElementType T;

    stringstream spy_filename;
    spy_filename << filename << ".gps";
    ofstream gps(spy_filename.str().c_str());
    stringstream content;

    int n = numRows(A.impl()), m = numCols(A.impl());
    DenseVector<Array<T> > e(m), y(m);

    // get first entry
    e(1) = 1.;
    y = A.impl()*e;
    T min, max;
    if (absoluteValues==true) {
        min = absolute(y(1)); max = absolute(y(1));
    } else {
        min = y(1); max = y(1);
    }
    e(1) = 0.;

    int nnz = 0;
    T value = 0.0;
    for (int j=1; j<=m; ++j) {
        e(j) = 1.;
        y = A.impl()*e;

        for (int i = 1; i<=n; ++i) {
            if (absoluteValues==true) {
                value = absolute(y(i));
            } else {
                value = y(i);
            }
            if (absolute(value)>=eps) {
                ++nnz;
                if (value>max) {
                    max = value;
                } else if (value<min) {
                    min = value;
                }
                // i=column => y-value, j=row => x-value
                content << j << " " << i << " " << value << endl;
            }
        }
        e(j) = 0.;
    }

    gps << "reset; set term png; set output '" << filename << ".png'" << endl;
    if (n==m) {
        gps << "set size square" << endl;
        gps << "set title '$n = " << n << ", n^2 = " << n*n;
    } else {
        gps << "set size ratio " << T(n)/T(m) << endl;
        gps << "set title '$n = " << n << ", m = " << m << ", n*m = " << n*m;
    }

    gps << ", nnz = " << nnz << "$'" << endl
        << "set xrange [1:" << n << "]" << endl
        << "set yrange [1:" << m << "] reverse" << endl
        << "unset xtics; unset ytics; unset ztics" << endl
        << "set pm3d map" << endl
        << "set cbrange[" << min << ":" << max << "]" << endl;

    gps << "set palette model RGB" << endl;
    if (absoluteValues==true) {
        T cbtics = exp(1.0/6.0 * log(max/min));
        gps << "set logscale cb" << endl;
        gps << "set cbtics " << min << ", " << cbtics << ", " << max << endl;
    } else {
        T cbtics = (max-min)/6.0;
        gps << "set cbtics " << min << ", " << cbtics << ", " << max << endl;
    }
    gps << "splot '-' with points pointtype 7 pointsize 0.5 palette notitle"
        << endl << content.str().c_str() << endl;
    gps << "set output" << endl
        << "set term pop" << endl;

    gps.close();
}
} // namespace lawa

