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

#include <list>

namespace lawa {

template <typename T, typename Basis>
void
getSingularPoints(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff, DenseVector<Array<T> > &sing_pts)
{

    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;
    std::list<T> temp;
    for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
        
        DenseVector<Array<T> > bf_singpts = basis.generator((*it).first.xtype).singularSupport((*it).first.j, (*it).first.k);

        for (int i = bf_singpts.firstIndex(); i <= bf_singpts.lastIndex(); ++i) {
            temp.push_back(bf_singpts(i));
        }
    }
    temp.sort(); temp.unique();
    sing_pts.engine().resize((int)temp.size());
    int i = 1;
    for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it ) {
        sing_pts(i) = *it; ++i;
    }
}

template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     const Preconditioner &P, T (*u)(T), const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());

    DenseVector<Array<T> > sing_pts;
    getSingularPoints(basis, coeff, sing_pts);

    for (int i=sing_pts.firstIndex(); i<=sing_pts.lastIndex(); ++i) {
        T x = sing_pts(i);
        T appr = 0.0;
        T exact= u(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j;
            long k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            
            appr += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);

        }
        plotfile << x << " " << exact << " " << appr  << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     Preconditioner &P, T (*u)(T), T (*du)(T), T a, T b, T h, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=0., d_appr = 0.0;
        T exact= u(x);
        T d_exact= du(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j;
            long k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            
            appr   += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);
            d_appr += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,1);

        }
        plotfile << x << " " << exact << " " << d_exact << " " << appr << " " << d_appr << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis, typename Preconditioner>
void
w_plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
       const Preconditioner &P, T (*u)(T), T (*w)(T), T a, T b, T h, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=0.;
        T exact= u(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j;
            long k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            appr   += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);
        }
        plotfile << x << " " << exact << " " << exact*w(x)
                      << " " << appr << " " << appr*w(x) << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       Preconditioner &P, T (*u)(T,T), T a1, T b1, T a2, T b2, T h, const char* filename)
{

    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(16);

    for (T x=a1; x<=b1; x+=h) {
        for (T y=a2; y<=b2; y+=h) {
            T appr = 0.0;
            T exact= u(x,y);
            for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype;
                XType xtype_y = (*it).first.index2.xtype;
                int j_x = (*it).first.index1.j, j_y = (*it).first.index2.j;
                long k_x = (*it).first.index1.k, k_y = (*it).first.index2.k;

                T coeff = (*it).second, prec = P((*it).first);
                
                appr += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,0);

            }
            plotfile << x << " " << y << " " << exact << " " << appr  << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       const Preconditioner &P, T (*u)(T,T), T (*dy_u)(T,T), T a1, T b1, T a2, T b2, 
       T h1, T h2, const char* filename)
{

    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());

    for (T x=a1; x<=b1; x+=h1) {
        for (T y=a2; y<=b2; y+=h2) {
            T appr = 0.0;
            T dy_appr = 0.0;
            T exact= u(x,y);
            T dy_exact = dy_u(x,y);
            for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype;
                XType xtype_y = (*it).first.index2.xtype;
                int j_x = (*it).first.index1.j, j_y = (*it).first.index2.j;
                long k_x = (*it).first.index1.k, k_y = (*it).first.index2.k;

                T coeff = (*it).second, prec = P((*it).first);
                
                appr    += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,0);
                dy_appr += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,1);

            }
            plotfile << x << " " << y << " " << exact << " " << appr  << " "
                     << dy_exact << " " << dy_appr << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis>
void
plotCoeff(const Coefficients<Lexicographical,T,Index1D > &coeff, const Basis &basis,
          const char* filename, bool locally_single_scale, bool interval)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_it;
    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    int shift = locally_single_scale ? 0 : 1;

    int j0  =  100000;
    int J   = -100000;
    T left  =  100000.;
    T right = -100000.;
    T maxCoeff = 0.;

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        int  j = (*it).first.j;
        long k = (*it).first.k;

        j0 = std::min(j0, j);
        J  = std::max(J, j);
        if ((*it).first.xtype == XBSpline) {
            maxCoeff = std::max(maxCoeff, (T)fabs((*it).second));
            if (locally_single_scale) {
                left  = std::min(left,  basis.refinementbasis.mra.phi.support(j,k).l1);
                right = std::max(right, basis.refinementbasis.mra.phi.support(j,k).l2);
            }
            else {
                left  = std::min(left,  basis.mra.phi.support(j,k).l1);
                right = std::max(right, basis.mra.phi.support(j,k).l2);
            }
        }
        else {
            maxCoeff = std::max(maxCoeff, (T)fabs((*it).second));
            left  = std::min(left, basis.psi.support(j,k).l1);
            right = std::max(right, basis.psi.support(j,k).l2);
        }
    }

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps enh color; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    int i=1;
    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        T lineWidth = 0.1;
        T ctr, fromX, toX, fromY, toY;
        T color = 0.0;
        int j       = (*it).first.j;
        long k  = (*it).first.k;
        XType xtype = (*it).first.xtype;

        if (xtype==XBSpline) {
            if (!interval) {
                fromX = basis.mra.phi.support(j,k).l1;
                toX   = basis.mra.phi.support(j,k).l2;
            }
            else {
                if (locally_single_scale) {
                    T h = 1./T(basis.refinementbasis.mra.cardI(j));
                    fromX  = std::max( ( k-basis.refinementbasis.mra.rangeI(j).firstIndex() ) * h, (T)0.);
                    toX    = std::min( ( k-basis.refinementbasis.mra.rangeI(j).firstIndex() + 1 ) * h, (T)1.);
                }
                else {
                    T h = 1./T(basis.mra.cardI(j));
                    fromX  = std::max( ( k-basis.mra.rangeI(j).firstIndex() ) * h, (T)0.);
                    toX    = std::min( ( k-basis.mra.rangeI(j).firstIndex() + 1 ) * h, (T)1.);
                }
            }
            fromY = j-shift-0.5;
            toY   = j-shift+0.5;
            color = fabs((*it).second) / maxCoeff;
        }

        else {
            if (!interval) {
                fromX = basis.psi.support(j,k).l1;
                toX   = basis.psi.support(j,k).l2;
            }
            else {
                T h = 1./T(basis.cardJ(j));
                fromX  = std::max( ( k-basis.rangeJ(j).firstIndex() ) * h, (T)0.);
                toX    = std::min( ( k-basis.rangeJ(j).firstIndex() + 1 ) * h, (T)1.);
            }
            fromY = j-0.5;
            toY   = j+0.5;
            color = fabs((*it).second) / maxCoeff;
        }
        gps << "set arrow " << i << " from " << fromX << ", " << fromY
                                 << " to " << fromX << ", " << toY << "lt -1 lw 3 nohead" << std::endl;
        ++i;
        gps << "set arrow " << i << " from " << fromX << ", " << fromY
                                 << " to " << toX << ", " << fromY << "lt -1 lw 3 nohead" << std::endl;
        ++i;
        gps << "set arrow " << i << " from " << toX << ", " << fromY
                                 << " to " << toX << ", " << toY << " lt -1 lw 3 nohead" << std::endl;
        ++i;
        gps << "set arrow " << i << " from " << fromX << ", " << toY
                                 << " to " << toX << ", " << toY << " lt -1 lw 3 nohead" << std::endl;
        ++i;
        /*
        gps << "set object rectangle from " << fromX << ", " << fromY
            << " to " << toX << "," << toY << " fc rgb";
        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";
        gps << " linewidth " << 0.1 << " fillstyle solid" << std::endl;
        */
    }


    gps << "set xrange["<< left <<":"<< right <<"]" << std::endl;
    gps << "set yrange[" << j0-shift-0.5 << ":" << J+0.5 << "]" << std::endl;
    //gps << "set xtics 1" << endl;
    //gps << "set ytics ('" << j0 << "' " << j0-1;
    //gps << ", '" << j0 << "' " << j0;
    //for (int j = j0+1; j <= J; ++j) {
    //    gps << ", '" << j << "' " << j;
    //}
    //gps << ")" << std::endl;
    gps << "plot " << j0-shift-0.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();
}

template <typename T, typename Basis>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff, const Basis &basis, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator const_it;
    if (coeff.size() == 0) {
        return;
    }
    std::cout << "plotCoeff was called!" << std::endl;

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps enh color; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    T maxCoeffSca = -1.0;
    T maxCoeffWav = -1.0;
    int j = (*coeff.begin()).second.j, k = (*coeff.begin()).second.k;
    int j0 = j;
    int J  = j;
    T a_sca = 5000.0, a_wav = 5000.0;
    T b_sca = -5000.0, b_wav = -5000.0;
    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        j = (*it).second.j; k = (*it).second.k;
        j0 = std::min(j0, j);
        J  = std::max(J, j);
        if ((*it).second.xtype == XBSpline) {
            maxCoeffSca = std::max(maxCoeffSca, fabs((*it).first));
            a_sca = std::min(a_sca, basis.mra.phi.support(j,k).l1);
            b_sca = std::max(b_sca, basis.mra.phi.support(j,k).l2);
        }
        else {
            maxCoeffWav = std::max(maxCoeffWav, fabs((*it).first));
            a_wav = std::min(a_wav, basis.psi.support(j,k).l1);
            b_wav = std::max(b_wav, basis.psi.support(j,k).l2);
        }
    }
    T maxCoeff = std::max(maxCoeffWav,maxCoeffSca);
    T l1_sca = basis.mra.phi.support(0,0).l1, l2_sca = basis.mra.phi.support(0,0).l2;
    T l1_wav = basis.psi.support(0,0).l1,     l2_wav = basis.psi.support(0,0).l2;

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        T lineWidth = 0.1;
        T ctr, fromX, toX, fromY, toY;
        T color = 0.0;

        if ((*it).second.xtype==XBSpline) {
            long k1 = ceil(pow2i<T>((*it).second.j)*a_sca - l1_sca), k2 = floor(pow2i<T>((*it).second.j)*b_sca - l2_sca);
            int N = k2 - k1 + 1;
            fromX = a_sca + ((*it).second.k-k1)*(b_sca-a_sca)/(T)N;
            toX   = a_sca + ((*it).second.k-k1+1)*(b_sca-a_sca)/(T)N;

            fromY = (*it).second.j-1.5;
            toY   = (*it).second.j-0.5;
            color = fabs((*it).first) / maxCoeffSca;
        }

        else {
            long k1 = ceil(pow2i<T>((*it).second.j)*a_wav - l1_wav), k2 = floor(pow2i<T>((*it).second.j)*b_wav - l2_wav);
            long N = k2 - k1 + 1;
            fromX = a_wav + ((*it).second.k-k1)*(b_wav-a_wav)/(T)N;
            toX   = fromX + std::max((b_wav-a_wav)/N,0.1);  //was 0.05

            fromY = (*it).second.j-0.5;
            toY   = (*it).second.j+0.5;
            color = fabs((*it).first) / maxCoeffWav;

        }

        gps << "set object rectangle from " << fromX << ", " << fromY
            << " to " << toX << "," << toY << " fc rgb";
        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";
        gps << " linewidth " << 0.1 << " fillstyle solid" << std::endl;

    }

    gps << "set xrange["<< std::min(a_sca, a_wav) <<":"<< std::max(b_sca, b_wav) <<"]" << std::endl;
    gps << "set yrange[" << j0-1.5 << ":" << J+0.5 << "]" << std::endl;
    //gps << "set xtics 1" << endl;
    gps << "set ytics ('" << j0 << "' " << j0-1;
    gps << ", '" << j0 << "' " << j0;
    for (int j = j0+1; j <= J; ++j) {
        gps << ", '" << j << "' " << j;
    }
    gps << ")" << std::endl;
    gps << "plot " << j0-1.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();
}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x, const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_coeff_abs_it;

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    T h = 5e-2;

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps color enh; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    const_coeff_abs_it first_element = coeff.begin();
    T max_value = fabs( (*first_element).first );
    T min_x=10000., max_x=-10000., min_y=10000., max_y=-10000.;
    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, j2=(*it).second.index2.j;
        long k1=(*it).second.index1.k, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;
        
        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        min_x = std::min(min_x,x); max_x = std::max(max_x,x);
        min_y = std::min(min_y,y); max_y = std::max(max_y,y);
    }

    T ratio = (max_y-min_y)/(max_x-min_x);

    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, k1=(*it).second.index1.k, j2=(*it).second.index2.j, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;
     
        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        min_x = std::min(min_x,x); max_x = std::max(max_x,x);
        min_y = std::min(min_y,y); max_y = std::max(max_y,y);
        T color = fabs((*it).first)/max_value;

        if (ratio<1) {
            gps << "set object rectangle from " << x-h << ", " << y-h/ratio << " to " << x+h << "," << y+h/ratio << "fc rgb ";
        }
        else {
            gps << "set object rectangle from " << x-h/ratio << ", " << y-h << " to " << x+h/ratio << "," << y+h << "fc rgb ";
        }

        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";

        gps << " fs solid 1.0" << std::endl;
    }
    gps << "set xlabel 'x'" << std::endl;
    gps << "set ylabel 'y'" << std::endl;
    gps << "set xrange[" << min_x-h << ":" << max_x+h << "]" << std::endl;
    gps << "set yrange[" << min_y-h << ":" << max_y+h << "]" << std::endl;
    gps << "plot " << std::min(min_x-h,min_y-h) << " w l lc rgb 'black' notitle " << std::endl;
    gps.close();
}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_coeff_abs_it;


    std::stringstream dataFilename;
    dataFilename << filename << ".dat";
    std::ofstream data(dataFilename.str().c_str());

    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, j2=(*it).second.index2.j;
        long k1=(*it).second.index1.k, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;

        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        if(fabs((*it).first) > 0){
          data << x << " " << y << " " << (*it).second << " " << 1. << std::endl;
        }
    }
    data.close();

}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<Lexicographical,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;


    std::stringstream dataFilename;
    dataFilename << filename << ".dat";
    std::ofstream data(dataFilename.str().c_str());

    for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).first.index1.j, j2=(*it).first.index2.j;
        long k1=(*it).first.index1.k, k2=(*it).first.index2.k;
        XType type1=(*it).first.index1.xtype, type2=(*it).first.index2.xtype;

        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);


        data << x << " " << y << " " << (*it).second << " " << -1. << std::endl;

    }
    data.close();

}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D_interval(const Coefficients<Lexicographical,T,Index> &coeff, const Basis_x &basis_x,
                            const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;


    std::stringstream dataFilename;
    dataFilename << filename << ".dat";
    std::ofstream data(dataFilename.str().c_str());

    int offset_x = 1-basis_x.mra.rangeI(basis_x.j0).firstIndex();
    int offset_y = 1-basis_y.mra.rangeI(basis_y.j0).firstIndex();
    for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).first.index1.j, j2=(*it).first.index2.j;
        long k1=(*it).first.index1.k, k2=(*it).first.index2.k;
        XType type1=(*it).first.index1.xtype, type2=(*it).first.index2.xtype;

        T x=0., y=0.;
        if (type1==XBSpline) {  x = T(k1+offset_x) / (basis_x.mra.cardI(j1)+1);      }
        else {                  x = T((k1-1)*2+1) / (basis_x.mra.cardI(j1+1)+1);  }
        if (type2==XBSpline) {  y = T(k2+offset_y) / (basis_y.mra.cardI(j2)+1);      }
        else {                  y = T((k2-1)*2+1) / (basis_y.mra.cardI(j2+1)+1);  }

        //center of the support
        if(fabs((*it).second) > 0){
          data << x << " " << y << " " << (*it).second << " " << -1. << std::endl;
        }
    }
    data.close();

}

template <typename T, typename Basis>
void
plotScatterCoeff(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis,
                 const char* filename, bool useSupportCenter,  T left_x1, T right_x1,
                 T left_x2, T right_x2)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    std::cerr << "plotScatterCoeff: #supp u = " << coeff.size() << std::endl;
    std::stringstream datafilename;
    datafilename << filename << ".dat";
    std::ofstream datafile(datafilename.str().c_str());

    T RightmLeft_x1 = right_x1-left_x1, RightmLeft_x2 = right_x2-left_x2;

    if (useSupportCenter) {
        for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int  j_x=(*it).first.index1.j, j_y=(*it).first.index2.j;
            long k_x=(*it).first.index1.k,  k_y=(*it).first.index2.k;
            XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype;

            Support<T> supp_x = basis.first.generator(xtype_x).support(j_x,k_x);
            supp_x.l1 *= RightmLeft_x1; supp_x.l2 *= RightmLeft_x1;
            supp_x.l1 += left_x1;       supp_x.l2 += left_x1;
            Support<T> supp_y = basis.second.generator(xtype_y).support(j_y,k_y);
            supp_y.l1 *= RightmLeft_x2; supp_y.l2 *= RightmLeft_x2;
            supp_y.l1 += left_x2;       supp_y.l2 += left_x2;

            T x=0., y=0.;
            x = (supp_x.l2 + supp_x.l1)/(T)2.;
            y = (supp_y.l2 + supp_y.l1)/(T)2.;

            //center of the support
            //if(fabs((*it).second) > 0){
              datafile << x << " " << y << " " << (*it).second << " " << 0. << std::endl;
            //}
        }
    }
    else {
        int offset_x = 1-basis.first.mra.rangeI(basis.first.j0).firstIndex();
        int offset_y = 1-basis.second.mra.rangeI(basis.second.j0).firstIndex();
        for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int  j_x=(*it).first.index1.j, j_y=(*it).first.index2.j;
            long k_x=(*it).first.index1.k,  k_y=(*it).first.index2.k;
            XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype;

            T x=0., y=0.;
            if (xtype_x==XBSpline) {  x = T(k_x+offset_x) / (basis.first.mra.cardI(j_x)+1);      }
            else {                    x = T((k_x-1)*2+1) / (basis.first.mra.cardI(j_x+1)+1);  }
            if (xtype_y==XBSpline) {  y = T(k_y+offset_y) / (basis.second.mra.cardI(j_y)+1);      }
            else {                    y = T((k_y-1)*2+1) / (basis.second.mra.cardI(j_y+1)+1);  }

            x *= RightmLeft_x1; x += left_x1;
            y *= RightmLeft_x1; y += left_x1;
            //center of the support
            //if(fabs((*it).second) > 0){
              datafile << x << " " << y << " " << (*it).second << " " << 0. << std::endl;
            //}
        }
    }
    datafile.close();
}

template <typename T, typename Basis>
void
plotScatterCoeff(const Coefficients<Lexicographical,T,Index3D> &coeff, const Basis &basis,
                 const char* filename, bool useSupportCenter)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff_it;

    std::stringstream datafilename;
    datafilename << filename << ".dat";
    std::ofstream datafile(datafilename.str().c_str());

    if (useSupportCenter) {
        for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int  j_x=(*it).first.index1.j, j_y=(*it).first.index2.j,  j_z=(*it).first.index3.j;
            long k_x=(*it).first.index1.k,  k_y=(*it).first.index2.k, k_z=(*it).first.index3.k;
            XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype,
                  xtype_z=(*it).first.index3.xtype;

            Support<T> supp_x = basis.first.generator(xtype_x).support(j_x,k_x);
            Support<T> supp_y = basis.first.generator(xtype_y).support(j_y,k_y);
            Support<T> supp_z = basis.first.generator(xtype_z).support(j_z,k_z);

            T x=0., y=0., z=0.;
            x = (supp_x.l2 + supp_x.l1)/(T)2.;
            y = (supp_y.l2 + supp_y.l1)/(T)2.;
            z = (supp_z.l2 + supp_z.l1)/(T)2.;
            //center of the support
            if(fabs((*it).second) > 0){
              datafile << x << " " << y << " " << z << " " << (*it).second << " " << -1. << std::endl;
            }
        }
    }
    else {
        int offset_x = 1-basis.first.mra.rangeI(basis.first.j0).firstIndex();
        int offset_y = 1-basis.second.mra.rangeI(basis.second.j0).firstIndex();
        int offset_z = 1-basis.third.mra.rangeI(basis.third.j0).firstIndex();
        for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int  j_x=(*it).first.index1.j, j_y=(*it).first.index2.j,  j_z=(*it).first.index3.j;
            long k_x=(*it).first.index1.k,  k_y=(*it).first.index2.k, k_z=(*it).first.index3.k;
            XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype,
                  xtype_z=(*it).first.index3.xtype;

            T x=0., y=0., z=0.;
            if (xtype_x==XBSpline) {  x = T(k_x+offset_x) / (basis.first.mra.cardI(j_x)+1);      }
            else {                    x = T((k_x-1)*2+1) / (basis.first.mra.cardI(j_x+1)+1);  }
            if (xtype_y==XBSpline) {  y = T(k_y+offset_y) / (basis.second.mra.cardI(j_y)+1);      }
            else {                    y = T((k_y-1)*2+1) / (basis.second.mra.cardI(j_y+1)+1);  }
            if (xtype_z==XBSpline) {  z = T(k_z+offset_z) / (basis.third.mra.cardI(j_z)+1);      }
            else {                    z = T((k_z-1)*2+1) / (basis.third.mra.cardI(j_z+1)+1);  }

            //center of the support
            if(fabs((*it).second) > 0){
              datafile << x << " " << y << " " << z << " " << (*it).second << " " << -1. << std::endl;
            }
        }
    }
    datafile.close();
}


}  // namespace lawa

