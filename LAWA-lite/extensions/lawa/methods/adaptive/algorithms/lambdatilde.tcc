namespace lawa {

/*
 * Realizations of lambdaTilde1d
 */

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis, 
                  int s_tilde, int jmin, int jmax=100, bool update=false)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;

    int j = lambda.j, k = lambda.k;
    int d = psi.d;
    IndexSet<Index1D> ret;
    Support<T> support_refbspline = phi.support(0,0);
    Support<T> support_refwavelet = psi.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << phi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            BSpline<T,Primal,R,CDF> phi_row(d);
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (overlap(supp, phi.support(j,k_row)) > 0) {
                    //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi_row.support(j,k_row) << " " << supp  << std::endl;
                    ret.insert(Index1D(j,k_row,XBSpline));
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                T Pow2i_Mjrow = pow2i<T>(-j_row);
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support phi_col = " << singpts;
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                            if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) >= 0 ))) {
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp  << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if (overlap(supp, supp_row) > 0)  {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
        else {
            Support<T> supp = psi.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Bsplines with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = floor( pow2i<T>(jmin)*supp.l1 - phi.support(0,0).l2)-1;
                int kMax =  ceil( pow2i<T>(jmin)*supp.l2 - phi.support(0,0).l1)+1;
                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    if (    (overlap(supp, phi.support(jmin,k_row)) > 0)
                         && (!(distance(supp,phi.singularSupport(jmin,k_row)) >= 0 )) )
                    {
                        //std::cout << "lambdaTilde: BSpline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " " << supp  << std::endl;
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
                T Pow2i_Mjrow = pow2i<T>(-j_row);
                if (j_row>=j+1) {
                    DenseVector<Array<T> > singsupp = psi.optim_singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                            if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if ((overlap(supp, supp_row) > 0) && (!(distance(psi.singularSupport(j_row,k_row),supp) >= 0 ))) {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
/*
    // We assume that indices corresponding to non-zero values are already calculated up to level jmax-1 > jmin.
    // Therefore, we only have to add wavelets on level jmax.
    else {
        int j = lambda.j, k = lambda.k;
        int d = lambda.d, d_= lambda.d_;

        if (lambda.xtype == XBSpline) {
            assert(j == jmin);
            BSpline<T,Primal,R> phi_col(d);
            Support<T> supp = phi_col.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " <<  phi_col.support(j,k) << " " << phi_col.singularSupport(j,k) << endl;

            if (jmax <= j+s_tilde) {
                // Inserting all indices corresponding to Wavelets with intersecting support using
                // a) local compactness  b) matrix compression  c) vanishing moments
                Wavelet<T,Primal,R> psi_row(d,d_);
                DenseVector<Array<T> > singpts = phi_col.singularSupport(j,k);
                for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
                    int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
                    int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
                            //cout << "LambdaTilde: kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(j_row,k_row) << endl;
                            ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                        }
                    }
                }
            }
        }

        else {
            assert(j >= jmin);
            Wavelet<T,Primal,R> psi_col(d,d_);
            Support<T> supp = psi_col.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            Wavelet<T,Primal,R> psi_row(d,d_);
            if (jmax <= j+s_tilde) {
                int level_diff = 3;
                if (d == 2) level_diff = 2;
                if (d == 3) level_diff = 4;
                if (jmax >= j+level_diff) {                                                // level difference has to be large enough for vanishing entries due to vanishing moments
                    DenseVector<Array<T> > singpts = psi_col.singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                    for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
                        int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
                        int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
                                //cout << "LambdaTilde update: jmax = " << jmax << ", kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(jmax,k_row) << endl;
                                ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = ceil( pow2i(jmax)*supp.l1 - psi_row.support(0,0).l2);
                    int kMax = floor(pow2i(jmax)*supp.l2 - psi_row.support(0,0).l1);
                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        if ((overlap(supp, psi_row.support(jmax,k_row)) == true) && !((distSingularSupport(psi_row,psi_col,jmax,k_row,j,k) > 0.0 ))){
                            //cout << "lambdaTilde_R: Wavelet (" << j_row << ", " << k_row << "): " << psi_row.support(j_row,k_row) << endl;
                            ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                        }
                    }
                }
            }
        }


    }
    */
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Orthogonal,R,Multi> &basis,
                  int s_tilde, int jmin, int jmax=100, bool update=false)
{
    IndexSet<Index1D> ret;
    if (update) return ret;

    const BSpline<T,Orthogonal,R,Multi> &phi = basis.mra.phi;
    const Wavelet<T,Orthogonal,R,Multi> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;

    int j  = lambda.j;
    long k = lambda.k;
    int d  = basis.d;

    Support<T> max_support_refbspline = phi.max_support();
    Support<T> max_support_refwavelet = psi.max_support();

    //std::cerr << "numSplines = " << numSplines << ", numWavelets = " << numWavelets
    //          << ", max_support_refbspline = " << max_support_refbspline
    //          << ", max_support_refwavelet = " << max_support_refwavelet << std::endl;


    if (lambda.xtype == XBSpline) {

        Support<T> supp = phi.support(j,k);

        //std::cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << supp << std::endl;

        // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
        for (long k_row=k-(d-1)*numSplines; k_row<=k+(d-1)*numSplines; ++k_row) {
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi.support(j,k_row) << std::endl;
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {
            T Pow2i_Mjrow = pow2i<T>(-j_row);

            if (j_row>=j+2) {
                DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    long kMin = floor(pow2i<long double>(j_row)*singsupp(i) - max_support_refwavelet.l2)-1;
                    long kMax =  ceil(pow2i<long double>(j_row)*singsupp(i) - max_support_refwavelet.l1)+1;
                    // todo: revise these loops and test!!
                    for (long k_help=kMin; k_help<=kMax; ++k_help) {
                        for (long k_row=(k_help-1)*numWavelets+1; k_row<=k_help*numWavelets; ++k_row) {
                            Support<T> supp_row = psi.support(j_row,k_row);
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                            if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) >= 0 ))) {
                                //Attention: cast from long to int here if old index class is used!!
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
            }
            else {

                long kMin = floor(pow2i<long double>(j_row)*supp.l1 - max_support_refwavelet.l2)-1;
                long kMax =  ceil(pow2i<long double>(j_row)*supp.l2 - max_support_refwavelet.l1)+1;
                // todo: revise these loops and test!!
                for (long k_help=kMin; k_help<=kMax; ++k_help) {
                    for (long k_row=(k_help-1)*numWavelets+1; k_row<=k_help*numWavelets; ++k_row) {
                        Support<T> supp_row = psi.support(j_row,k_row);
                        //std::cout << "LambdaTilde: Wavelet (" << j_row << ", " << k_row << "): " << supp_row << " " << supp << ", overlap=" << overlap(supp, psi.support(j_row,k_row)) << std::endl;
                        if (overlap(supp, supp_row) > 0)  {
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        //std::cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << supp << std::endl;
        //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

        // Inserting all indices corresponding to Bsplines with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        if (fabs(j - jmin) <= s_tilde) {
            long kMin = floor( pow2i<long double>(jmin)*supp.l1 - max_support_refbspline.l2);
            long kMax =  ceil( pow2i<long double>(jmin)*supp.l2 - max_support_refbspline.l1);
            // todo: revise these loops and test!!
            for (long k_help=kMin; k_help<=kMax; ++k_help) {
                for (long k_row=(k_help-1)*numSplines+1; k_row<=k_help*numSplines; ++k_row) {
                    if (    (overlap(supp, phi.support(jmin,k_row)) > 0)
                         && (!(distance(supp,phi.singularSupport(jmin,k_row)) >= 0 )) )
                    {
                        //std::cout << "lambdaTilde: BSpline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " " << distance(supp,phi.singularSupport(jmin,k_row))  << std::endl;
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                }
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
            long double Pow2i_Mjrow = pow2i<long double>(-j_row);
            if (j_row>=j+4) {
                DenseVector<Array<T> > singsupp = psi.singularSupport(j,k);

                //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    long kMin = floor(pow2i<long double>(j_row)*singsupp(i) - max_support_refwavelet.l2);
                    long kMax =  ceil(pow2i<long double>(j_row)*singsupp(i) - max_support_refwavelet.l1);
                    // todo: revise these loops and test!!
                    for (long k_help=kMin; k_help<=kMax; ++k_help) {
                        for (long k_row=(k_help-1)*numWavelets+1; k_row<=k_help*numWavelets; ++k_row) {
                            Support<T> supp_row = psi.support(j_row,k_row);
                            if (   (overlap(supp, supp_row) > 0)
                                && (!(distance(singsupp,supp_row) >= 0 ))){
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
            }
            else {
                long kMin = floor(pow2i<long double>(j_row)*supp.l1 - max_support_refwavelet.l2)-1;
                long kMax =  ceil(pow2i<long double>(j_row)*supp.l2 - max_support_refwavelet.l1)+1;
                // todo: revise these loops and test!!
                for (long k_help=kMin; k_help<=kMax; ++k_help) {
                    for (long k_row=(k_help-1)*numWavelets+1; k_row<=k_help*numWavelets; ++k_row) {
                        Support<T> supp_row = psi.support(j_row,k_row);
                        //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << std::endl;
                        if ((overlap(supp, supp_row) > 0)
                            && (!(distance(psi.singularSupport(j_row,k_row),supp) >= 0 ))) {

                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Orthogonal,Interval,Multi> &basis,
                  int s_tilde, int jmin, int jmax=100, bool update=false)
{
    IndexSet<Index1D> ret;
    if (update) return ret;

    if (basis.j0!=jmin) {
        std::cerr << "lambdaTilde1d_PDE: Something is wrong with the minimal level: "
                  << basis.j0 << "!=" << jmin << std::endl;
        exit(1);
    }

    const BSpline<T,Orthogonal,Interval,Multi> &phi = basis.mra.phi;
    const Wavelet<T,Orthogonal,Interval,Multi> &psi = basis.psi;
    Support<T> max_support(-1.,1.); //maximum support of a generator in a multiwavelet basis.
    // Right boundary has less associated boundary functions than left boundary.
    int numWavelets = (int)std::max(basis._numLeftParts,basis._numInnerParts);

    int j = lambda.j;
    long k = lambda.k;
    //int d = psi.d;


    if (lambda.xtype == XBSpline) {

        Support<T> supp = phi.support(j,k);

        //std::cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << supp << std::endl;

        // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
        for (long k_row=basis.mra.rangeI(j).firstIndex(); k_row<=basis.mra.rangeI(j).lastIndex(); ++k_row) {
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi.support(j,k_row) << std::endl;
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {
            T Pow2i_Mjrow = pow2i<T>(-j_row);

            // Check if boundary wavelets (not all of them have vanishing moments) intersect
            for (int k_row=basis.rangeJL(j_row).firstIndex(); k_row<=basis.rangeJL(j_row).lastIndex(); ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                if (overlap(supp_row, supp) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
            for (int k_row=basis.rangeJR(j_row).firstIndex(); k_row<=basis.rangeJR(j_row).lastIndex(); ++k_row) {
                if (j_row < 0) {
                    std::cerr << "ERROR: no negative levels allowed here!" << std::endl;
                    exit(1);
                }
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                if (overlap(supp_row, supp) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }


            if (j_row>=j+4) {
                DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    long kMin = floor(pow2i<long>(j_row)*singsupp(i) - max_support.l2)-1;
                    long kMax =  ceil(pow2i<long>(j_row)*singsupp(i) - max_support.l1)+1;
                    long k_first = std::max((kMin-1)*numWavelets+1, (long)basis.rangeJI(j_row).firstIndex());
                    long k_last  = std::min(kMax*numWavelets,       (long)basis.rangeJI(j_row).lastIndex());
                    for (long k_row=k_first; k_row<=k_last; ++k_row) {
                        Support<T> supp_row = psi.support(j_row,k_row);
                        //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                        if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) >= 0 ))) {
                            //Attention: cast from long to int here if old index class is used!!
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
            else {
                for (int j_row=j; j_row<=j+3; ++j_row) {
                    for (int k_row=basis.rangeJI(j_row).firstIndex(); k_row<=basis.rangeJI(j_row).lastIndex(); ++k_row) {
                        Support<T> supp_row = psi.support(j_row,k_row);
                        if (overlap(supp, supp_row) > 0)  {
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        bool hasVanishingMoments = ((k>basis.rangeJ(j).lastIndex()) && (k<basis.rangeJ(j).firstIndex())) ? true : false;
        // Inserting all indices corresponding to Bsplines with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        if (fabs(j - jmin) <= s_tilde) {
            for (long k_row=basis.mra.rangeI(jmin).firstIndex(); k_row<=basis.mra.rangeI(jmin).lastIndex(); ++k_row) {
                if (overlap(supp, phi.support(jmin,k_row)) > 0) {
                    if (hasVanishingMoments) {
                        if (!(distance(phi.singularSupport(jmin,k_row),supp) >= 0 )) {
                            ret.insert(Index1D(jmin,k_row,XBSpline));
                        }
                    }
                    else {
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                }
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {

            // Check if boundary wavelets (not all of them have vanishing moments) intersect
            for (int k_row=basis.rangeJL(j_row).firstIndex(); k_row<=basis.rangeJL(j_row).lastIndex(); ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                if (overlap(supp_row, supp) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
            for (int k_row=basis.rangeJR(j_row).firstIndex(); k_row<=basis.rangeJR(j_row).lastIndex(); ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << supp_row << std::endl;
                if (overlap(supp_row, supp) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }

            if (j_row>=j+4) {
                DenseVector<Array<T> > singsupp = psi.singularSupport(j,k);

                //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    long kMin = floor(pow2i<long double>(j_row)*singsupp(i) - max_support.l2);
                    long kMax =  ceil(pow2i<long double>(j_row)*singsupp(i) - max_support.l1);
                    long k_first = std::max((kMin-1)*numWavelets,(long)basis.rangeJI(j_row).firstIndex());
                    long k_last  = std::min(kMax*numWavelets,    (long)basis.rangeJI(j_row).lastIndex());

                    for (long k_row=k_first; k_row<=k_last; ++k_row) {
                        Support<T> supp_row = psi.support(j_row,k_row);
                        if (   (overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
            else {
                long kMin = floor(pow2i<long double>(j_row)*supp.l1 - max_support.l2)-1;
                long kMax =  ceil(pow2i<long double>(j_row)*supp.l2 - max_support.l1)+1;
                long k_first = std::max((kMin-1)*numWavelets,(long)basis.rangeJI(j_row).firstIndex());
                long k_last  = std::min(kMax*numWavelets,    (long)basis.rangeJI(j_row).lastIndex());

                for (long k_row=k_first; k_row<=k_last; ++k_row) {
                    Support<T> supp_row = psi.support(j_row,k_row);
                    //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << std::endl;
                    if (overlap(supp, supp_row) > 0) {
                        if (hasVanishingMoments) {
                            if (!(distance(psi.singularSupport(j_row,k_row),supp) >= 0 )) {
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                        else {
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,R,SparseMulti> &basis,
                  int /*s_tilde*/, int jmin, int jmax=100, bool update=false)
{
    IndexSet<Index1D> ret;
    if (update) return ret;

    jmax = std::min(jmax,JMAX);

    const BSpline<T,Primal,R,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,R,SparseMulti> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;

    int  j = lambda.j;
    long k = lambda.k;
    int  d = psi.d;
    assert(d==4);

    Support<T> max_support_refbspline = phi.max_support();
    Support<T> max_support_refwavelet = psi.max_support();

    if (lambda.xtype == XBSpline) {

        Support<T> supp = phi.support(j,k);

        long kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refbspline.l2) - 1;
        long kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refbspline.l1) + 1;

        for (long k_row=kMin*numSplines; k_row<=kMax*numSplines; ++k_row) {
            //std::cerr << "  -> bspline (" << j << ", " << k_row << "): " << phi.support(j,k_row) << " vs. " << supp << std::endl;
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }

        kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;
        kMin = kMin*numWavelets;
        kMax = kMax*numWavelets;


        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j,k_row);
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
            if (overlap(supp, supp_row) > 0)  {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        if (j==jmin) {

            long kMin = floor( pow2i<long double>(jmin)*supp.l1 - max_support_refbspline.l2) - 1;
            long kMax =  ceil( pow2i<long double>(jmin)*supp.l2 - max_support_refbspline.l1) + 1;

            kMin = kMin*numSplines;
            kMax = kMax*numSplines;

            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                //std::cerr << "  -> bspline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " vs. " << supp << std::endl;
                if (overlap(supp, phi.support(jmin,k_row)) > 0) {
                    ret.insert(Index1D(jmin,k_row,XBSpline));
                }
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=std::max(j-1,jmin); j_row<=std::min(j+1,jmax); ++j_row) {

            long kMin = floor( pow2i<long double>(j_row)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
            long kMax =  ceil( pow2i<long double>(j_row)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;

            kMin = kMin*numWavelets;
            kMax = kMax*numWavelets;


            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cerr << "  -> wavelet (" << j_row << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
                if (overlap(supp, supp_row) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,RPlus,SparseMulti> &basis,
                  int /*s_tilde*/, int jmin, int jmax=100, bool update=false)
{
    IndexSet<Index1D> ret;
    if (update) return ret;

    jmax = std::min(jmax,JMAX);

    const BSpline<T,Primal,RPlus,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,RPlus,SparseMulti> &psi = basis.psi;
    int numSplines = (int)phi._numSplines;
    int numWavelets = (int)psi._numSplines;

    int  j = lambda.j;
    long k = lambda.k;
    int  d = psi.d;
    assert(d==4);

    Support<T> max_support_refbspline = phi.max_support();
    Support<T> max_support_refwavelet = psi.max_support();

    if (lambda.xtype == XBSpline) {
        Support<T> supp = phi.support(j,k);

        long kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refbspline.l2) - 1;
        long kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refbspline.l1) + 1;
        kMin = std::max(kMin*numSplines, basis.mra.rangeIL(j).firstIndex());
        kMax = kMax*numSplines;
        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            //std::cerr << "  -> bspline (" << j << ", " << k_row << "): " << phi.support(j,k_row) << " vs. " << supp << std::endl;
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }

        kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;
        kMin = std::max(kMin*numWavelets, basis.rangeJL(j).firstIndex());
        kMax = kMax*numWavelets;
        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j,k_row);
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
            if (overlap(supp, supp_row) > 0)  {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        if (j==jmin) {
            long kMin = floor( pow2i<long double>(jmin)*supp.l1 - max_support_refbspline.l2) - 1;
            long kMax =  ceil( pow2i<long double>(jmin)*supp.l2 - max_support_refbspline.l1) + 1;
            kMin = std::max(kMin*numSplines, basis.mra.rangeIL(j).firstIndex());
            kMax = kMax*numSplines;
            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                if (overlap(supp, phi.support(jmin,k_row)) > 0) {
                    ret.insert(Index1D(jmin,k_row,XBSpline));
                }
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=std::max(j-1,jmin); j_row<=std::min(j+1,jmax); ++j_row) {
            long kMin = floor( pow2i<long double>(j_row)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
            long kMax =  ceil( pow2i<long double>(j_row)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;
            kMin = std::max(kMin*numWavelets, basis.rangeJL(j_row).firstIndex());
            kMax = kMax*numWavelets;
            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cerr << "  -> wavelet (" << j_row << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
                if (overlap(supp, supp_row) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Interval,SparseMulti> &basis,
                  int /*s_tilde*/, int jmin, int jmax=100, bool update=false)
{
    IndexSet<Index1D> ret;
    if (update) return ret;

    jmax = std::min(jmax,JMAX);

    const BSpline<T,Primal,Interval,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,Interval,SparseMulti> &psi = basis.psi;
    int numSplines = (int)basis.mra._numSplines;
    int numWavelets = (int)basis._numSplines;

    int  j = lambda.j;
    long k = lambda.k;
    int  d = psi.d;

    Support<T> max_support_refbspline = basis.mra.max_support();
    Support<T> max_support_refwavelet = basis.max_support();

    if (lambda.xtype == XBSpline) {
        Support<T> supp = phi.support(j,k);

        long kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refbspline.l2) - 1;
        long kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refbspline.l1) + 1;
        kMin = std::max(kMin*numSplines, basis.mra.long_rangeI(j).firstIndex());
        kMax = std::min(kMax*numSplines, basis.mra.long_rangeI(j).lastIndex());
        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            //std::cerr << "  -> bspline (" << j << ", " << k_row << "): " << phi.support(j,k_row) << " vs. " << supp << std::endl;
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }

        kMin = floor( pow2i<long double>(j)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
        kMax =  ceil( pow2i<long double>(j)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;
        kMin = std::max(kMin*numWavelets, basis.long_rangeJ(j).firstIndex());
        kMax = std::min(kMax*numWavelets, basis.long_rangeJ(j).lastIndex());
        for (long k_row=kMin; k_row<=kMax; ++k_row) {
            Support<T> supp_row = psi.support(j,k_row);
            //std::cerr << "  -> wavelet (" << j << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
            if (overlap(supp, supp_row) > 0)  {
                ret.insert(Index1D(j,k_row,XWavelet));
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        if (j==jmin) {
            long kMin = floor( pow2i<long double>(jmin)*supp.l1 - max_support_refbspline.l2) - 1;
            long kMax =  ceil( pow2i<long double>(jmin)*supp.l2 - max_support_refbspline.l1) + 1;
            kMin = std::max(kMin*numSplines, basis.mra.long_rangeI(j).firstIndex());
            kMax = std::min(kMax*numSplines, basis.mra.long_rangeI(j).lastIndex());
            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                if (overlap(supp, phi.support(jmin,k_row)) > 0) {
                    ret.insert(Index1D(jmin,k_row,XBSpline));
                }
            }
        }
        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        for (int j_row=std::max(j-1,jmin); j_row<=std::min(j+1,jmax); ++j_row) {
            long kMin = floor( pow2i<long double>(j_row)*supp.l1 - max_support_refwavelet.l2) / 2 - 1;
            long kMax =  ceil( pow2i<long double>(j_row)*supp.l2 - max_support_refwavelet.l1) / 2 + 1;
            kMin = std::max(kMin*numWavelets, basis.long_rangeJ(j_row).firstIndex());
            kMax = std::min(kMax*numWavelets,  basis.long_rangeJ(j_row).lastIndex());
            for (long k_row=kMin; k_row<=kMax; ++k_row) {
                Support<T> supp_row = psi.support(j_row,k_row);
                //std::cerr << "  -> wavelet (" << j_row << ", " << k_row << "): " << supp_row << " vs. " << supp << std::endl;
                if (overlap(supp, supp_row) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Periodic,CDF> &basis, 
                  int s_tilde, int jmin, int jmax, bool update)
{
    BSpline<T,Primal,Periodic,CDF> phi_col(basis.mra), phi_row(basis.mra);
    Wavelet<T,Primal,Periodic,CDF> psi_col(basis), psi_row(basis);
    int j = lambda.j, k = lambda.k;

    IndexSet<Index1D> ret;
    Support<T> support_refbspline = phi_col.phiR.support(0,0);
    Support<T> support_refwavelet = psi_col.psiR.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi_col.phiR.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                int k_row_per = k_row;
                if(k_row_per < basis.mra.rangeI(j).firstIndex()){
                    k_row_per = basis.mra.rangeI(j).lastIndex() + ((1 - (basis.mra.rangeI(j).firstIndex() - k_row_per))%basis.mra.cardI(j));
                }
                if(k_row_per > basis.mra.rangeI(j).lastIndex()){
                    k_row_per = basis.mra.rangeI(j).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j));
                }
                if (overlap(phi_row.support(j,k_row_per), phi_col.support(j,k)) > 0) {
                    ret.insert(Index1D(jmin,k_row_per,XBSpline));
               }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi_col.phiR.singularSupport(j,k);
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            int k_row_per = k_row;
                            if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                                k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                               }
                               if(k_row_per> basis.rangeJ(j_row).lastIndex()){
                                   k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                               }
                               Support<T> supp_row = psi_row.support(j_row,k_row_per);
                            if (((overlap(supp_row, phi_col.support(j,k)) > 0)) &&
                                (!(distance(phi_col.singularSupport(j,k),supp_row) >= 0 ))) {
                                ret.insert(Index1D(j_row,k_row_per,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        int k_row_per = k_row;
                        if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                        }
                        if(k_row_per > basis.rangeJ(j_row).lastIndex()){
                            k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                        }
                        if (overlap(psi_row.support(j_row,k_row_per), phi_col.support(j,k)) > 0) {
                            ret.insert(Index1D(j_row,k_row_per,XWavelet));
                        }
                    }
                }
            }
        }

        else {
            Support<T> supp = psi_col.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2);
                int kMax =  ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1);
                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    int k_row_per = k_row;
                    if(k_row_per < basis.mra.rangeI(jmin).firstIndex()){
                        k_row_per = basis.mra.rangeI(jmin).lastIndex() + ((1 - (basis.mra.rangeI(jmin).firstIndex() - k_row_per))%basis.mra.cardI(jmin));
                    }
                    if(k_row_per > basis.mra.rangeI(jmin).lastIndex()){
                        k_row_per = basis.mra.rangeI(jmin).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(jmin).lastIndex()))%basis.mra.cardI(jmin));
                    }

                    if (overlap(phi_row.support(jmin,k_row_per), psi_col.support(j,k)) > 0) {
                        ret.insert(Index1D(jmin,k_row_per,XBSpline));
                    }
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = psi_col.singularSupport(j,k);
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            int k_row_per = k_row;
                            if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                                k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                            }
                            if(k_row_per> basis.rangeJ(j_row).lastIndex()){
                                k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                            }
                            Support<T> supp_row = psi_row.support(j_row,k_row_per);
                            if (((overlap(supp_row, psi_col.support(j,k)) > 0)) &&
                                (!(distance(psi_col.singularSupport(j,k),supp_row) >= 0 )) ) {
                                    ret.insert(Index1D(j_row,k_row_per,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        int k_row_per = k_row;
                        if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                        }
                        if(k_row_per > basis.rangeJ(j_row).lastIndex()){
                            k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                        }
                        if ((overlap(psi_row.support(j_row,k_row_per), psi_col.support(j,k)) > 0) &&
                            !(distance(psi_row.singularSupport(j_row,k_row_per),supp) >= 0 ) ) {
                           ret.insert(Index1D(j_row,k_row_per,XWavelet));
                        }
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T, Construction Cons>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis,
                  int s_tilde, int jmin, int jmax, bool /*update*/) 
{
        using std::min;
        using std::max;

        BSpline<T,Primal,Interval,Cons> phi_col(basis.mra), phi_row(basis.mra);
        Wavelet<T,Primal,Interval,Cons> psi_col(basis), psi_row(basis);
        int j = lambda.j, k = lambda.k;
        IndexSet<Index1D> ret;

        if (lambda.xtype==XBSpline) {

            Support<T> supp_col = phi_col.support(j,k);

            //Adding B-Splines
            int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
            int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(jmin)),kMin), kMax);
            //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
            while ((kStart-1 >= kMin) && (overlap(supp_col, phi_row.support(jmin,max(kStart-1, kMin)))>0)) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)),kMax), kMin);
            assert((overlap(supp_col, phi_col.support(jmin,kEnd))>0));
            while ((kEnd+1 <= kMax) && (overlap(supp_col, phi_row.support(jmin,min(kEnd+1,kMax)))>0)) {
                ++kEnd;
            }

            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                ret.insert(Index1D(jmin,k_row,XBSpline));
            }

            //Adding Wavelets
            for (int j_row=jmin; j_row<=min(jmin+s_tilde, jmax); ++j_row) {

                int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
                int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
                //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
                //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }
                for (int k_row=kStart; k_row<=kEnd; k_row++) {
                    Range<int> rangeL = basis.rangeJL(j_row);
                    Range<int> rangeR = basis.rangeJR(j_row);
                    if ( (k_row <= rangeL.lastIndex()) || (k_row >= rangeR.firstIndex()) ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                        continue;
                    }
                    if  (!(distance(phi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) {    //singsupp!
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
        else {
            Support<T> supp_col = psi_col.support(j,k);

            //Adding B-Splines
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
                int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(jmin)), kMin), kMax);
                //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,phi_row.support(jmin,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)), kMax), kMin);
                //assert((overlap(supp_col, phi_row.support(jmin,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,phi_row.support(jmin,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }

                for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                    if  (distance(psi_col.singularSupport(j,k),phi_row.support(jmin,k_row)) < 0 ) {        //singsupp!
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                    else {
                        if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                             (k >= basis.rangeJR(j).firstIndex() )     ) {
                                ret.insert(Index1D(jmin,k_row,XBSpline));
                        }
                    }
                }
            }

            //Adding Wavelets
            for (int j_row=max(j-s_tilde,jmin); j_row<=min(j+s_tilde,jmax); ++j_row) {

                int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
                int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
                //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
                //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }
                for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                    if (distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) < 0 ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                        continue;
                    }
                    else {
                        if ( (k_row <= basis.rangeJL(j_row).lastIndex() )  ||
                             (k_row >= basis.rangeJR(j_row).firstIndex() )     ) {
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                        continue;
                    }

                    if (distance(psi_row.singularSupport(j_row,k_row),supp_col) < 0 ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                    else {
                        if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                             (k >= basis.rangeJR(j).firstIndex() )     ) {
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }

                    /*
                    if  ( (!(distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) &&    //singsupp!
                          (!(distance(psi_row.singularSupport(j_row,k_row), supp_col) >= 0 )) ) {            //singsupp!
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                    */
                }
            }
        }
        return ret;
}


template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE_WO_XBSpline(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis, 
                              int s_tilde, int jmin, int jmax)
{
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;
    int j = lambda.j, k = lambda.k;
    IndexSet<Index1D> ret;
    Support<T> support_refwavelet = psi.support(0,0);
    Support<T> supp = psi.support(j,k);

    // Inserting all indices corresponding to Wavelets with intersecting support using
    // a) local compactness  b) matrix compression  c) vanishing moments
    for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
        T Pow2i_Mjrow = pow2i<T>(-j_row);
        if (j_row>=j+2) {
            DenseVector<Array<T> > singsupp = psi.optim_singularSupport(j,k);
            for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));
                    if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
        else {
            int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
            int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (T(Pow2i_Mjrow*(support_refwavelet.l1+k_row)) >= T(Pow2i_Mjrow*(support_refwavelet.l2+k_row))) {
                    std::cout << "Warning3: Translation indices too large!!" << std::endl;
                    continue;
                }
                Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));
                if ((overlap(supp, supp_row) > 0) && (!(distance(psi.optim_singularSupport(j_row,k_row),supp) > 0 ))) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE_WO_XBSpline(const Index1D &/*lambda*/, const Basis<T,Orthogonal,R,Multi> &/*basis*/,
                              int /*s_tilde*/, int /*jmin*/, int /*jmax*/)
{
    assert(0);
    std::cerr << "lambdaTilde1d_PDE_WO_XBSpline not implemented for "
              << "Basis<T,Orthogonal,R,Multi>." << std::endl;
    exit(1);
    IndexSet<Index1D> ret;
    return ret;
}


template <typename T, Construction Cons>
IndexSet<Index1D>
lambdaTilde1d_WeightedPDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis,
                          int s_tilde_level, int jmin, int jmax, int s_tilde_singsupp)
{
    using std::min;
    using std::max;

    if (s_tilde_singsupp==-1) {
        s_tilde_singsupp = (int)((basis.d-1.5)*s_tilde_level/(1.+basis.d_))+1;
    }

    BSpline<T,Primal,Interval,Cons> phi_col(basis.mra), phi_row(basis.mra);
    Wavelet<T,Primal,Interval,Cons> psi_col(basis), psi_row(basis);
    int j = lambda.j, k = lambda.k;
    IndexSet<Index1D> ret;

    if (lambda.xtype==XBSpline) {

        Support<T> supp_col = phi_col.support(j,k);

        //Adding B-Splines
        int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
        int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(jmin)),kMin), kMax);
        //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
        while ((kStart-1 >= kMin) && (overlap(supp_col, phi_row.support(jmin,max(kStart-1, kMin)))>0)) {
            --kStart;
        }
        int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)),kMax), kMin);
        assert((overlap(supp_col, phi_col.support(jmin,kEnd))>0));
        while ((kEnd+1 <= kMax) && (overlap(supp_col, phi_row.support(jmin,min(kEnd+1,kMax)))>0)) {
            ++kEnd;
        }

        for (int k_row=kStart; k_row<=kEnd; ++k_row) {
            ret.insert(Index1D(jmin,k_row,XBSpline));
        }

        //Adding Wavelets
        for (int j_row=jmin; j_row<=min(jmin+s_tilde_level, jmax); ++j_row) {

            int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
            int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
            //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
            while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
            //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
            while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                ++kEnd;
            }
            for (int k_row=kStart; k_row<=kEnd; k_row++) {
                Range<int> rangeL = basis.rangeJL(j_row);
                Range<int> rangeR = basis.rangeJR(j_row);
                if ( (k_row <= rangeL.lastIndex()) || (k_row >= rangeR.firstIndex()) ) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                    continue;
                }
                // Cannot use vanishing moments -> test only for disjunct supports
                if  (overlap(phi_col.support(j,k),psi_row.support(j_row,k_row)) > 0) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    else {
        Support<T> supp_col = psi_col.support(j,k);

        //Adding B-Splines
        if (fabs(j - jmin) <= s_tilde_level) {
            int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
            int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(jmin)), kMin), kMax);
            //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
            while (kStart-1>=kMin && overlap(supp_col,phi_row.support(jmin,max(kStart-1,kMin)))>0) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)), kMax), kMin);
            //assert((overlap(supp_col, phi_row.support(jmin,kEnd))>0));
            while (kEnd+1<=kMax && overlap(supp_col,phi_row.support(jmin,min(kEnd+1,kMax)))>0) {
                ++kEnd;
            }

            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                // Cannot use vanishing moments -> test only for disjunct supports
                if  (overlap(psi_col.support(j,k),phi_row.support(jmin,k_row)) > 0 ) {
                    ret.insert(Index1D(jmin,k_row,XBSpline));
                }
                else {
                    if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                         (k >= basis.rangeJR(j).firstIndex() )     ) {
                            ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                }
            }
        }

        //Adding Wavelets
        for (int j_row=max(j-s_tilde_level,jmin); j_row<=min(j+s_tilde_level,jmax); ++j_row) {

            int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
            int kStart = min(max(iceil<int>(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
            //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
            while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
            //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
            while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                ++kEnd;
            }
            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                // Cannot use vanishing moments -> test only for disjunct supports
                if (overlap(psi_col.support(j,k),psi_row.support(j_row,k_row)) > 0 ) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                    continue;
                }
                else {
                    if ( (k_row <= basis.rangeJL(j_row).lastIndex() )  ||
                         (k_row >= basis.rangeJR(j_row).firstIndex() )     ) {
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                    continue;
                }

                // Cannot use vanishing moments -> test only for disjunct supports
                if (overlap(psi_row.support(j_row,k_row),supp_col) > 0 ) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
                else {
                    if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                         (k >= basis.rangeJR(j).firstIndex() )     ) {
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
    }
    return ret;
      
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_WeightedPDE(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis,
                          int s_tilde_level, int jmin, int jmax, int s_tilde_singsupp)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;

    if (s_tilde_singsupp==-1) {
        s_tilde_singsupp = (int)((basis.d-1.5)*s_tilde_level/(1.+basis.d_))+4;
    }

    int j = lambda.j, k = lambda.k;
    //int d = psi.d;
    IndexSet<Index1D> ret;
    Support<T> support_refbspline = phi.support(0,0);
    Support<T> support_refwavelet = psi.support(0,0);

    if (lambda.xtype == XBSpline) {
        Support<T> supp = phi.support(j,k);
        //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << phi_col.singularSupport(j,k) << endl;

        // Inserting all indices corresponding to B-Splines with intersecting support using local compactness
        int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
        int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
        for (int k_row=kMin; k_row<=kMax; ++k_row) {
            if (overlap(supp, phi.support(j,k_row)) > 0) {
                //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi_row.support(j,k_row) << " " << supp  << std::endl;
                ret.insert(Index1D(j,k_row,XBSpline));
            }
        }

        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments for level diff>=s_tilde_singsupp
        for (int j_row=j; j_row<=std::min(j+s_tilde_level, jmax); ++j_row) {        // realization of matrix compression via level threshold
            T Pow2i_Mjrow = pow2i<T>(-j_row);
            if (j_row>j+s_tilde_singsupp) {
                DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                //cout << "LambdaTilde: Singular support phi_col = " << singpts;
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) >= 0 ))) {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp  << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
            else {
                int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                    if (overlap(supp, supp_row) > 0)  {
                        //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
    }
    else {
        Support<T> supp = psi.support(j,k);
        //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

        // Inserting all indices corresponding to Bsplines with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments
        if (fabs(j - jmin) <= s_tilde_level) {
            int kMin = floor( pow2i<T>(jmin)*supp.l1 - phi.support(0,0).l2)-1;
            int kMax =  ceil( pow2i<T>(jmin)*supp.l2 - phi.support(0,0).l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (    (overlap(supp, phi.support(jmin,k_row)) > 0)
                     &&  (!(distance(supp,phi.singularSupport(jmin,k_row)) >= 0 )) )  {
                    //std::cout << "lambdaTilde: BSpline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " " << supp  << std::endl;
                    ret.insert(Index1D(jmin,k_row,XBSpline));
                }
            }
        }

        // Inserting all indices corresponding to Wavelets with intersecting support using
        // a) local compactness  b) matrix compression  c) vanishing moments for level diff >= s_tilde_singsupp
        for (int j_row=std::max(j-s_tilde_level,jmin); j_row<=std::min(j+s_tilde_level,jmax); ++j_row) {
            T Pow2i_Mjrow = pow2i<T>(-j_row);
            if (j_row>j+s_tilde_singsupp) {
                DenseVector<Array<T> > singsupp = psi.optim_singularSupport(j,k);
                //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                    int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
            else {
                int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                    if ((overlap(supp, supp_row) > 0) && (!(distance(psi.optim_singularSupport(j_row,k_row),supp) >= 0 ))) {
                        //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
    }
    return ret;
}


} // namespace lawa

