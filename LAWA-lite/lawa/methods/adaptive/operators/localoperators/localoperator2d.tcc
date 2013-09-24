namespace lawa {

template <typename LocalOperator1, typename LocalOperator2>
LocalOperator2D<LocalOperator1, LocalOperator2>
::LocalOperator2D(LocalOperator1 &_localoperator1, LocalOperator2 &_localoperator2)
: localoperator1(_localoperator1), localoperator2(_localoperator2),
  trialBasis_x1(_localoperator1.trialBasis), testBasis_x1(_localoperator1.testBasis),
  trialBasis_x2(_localoperator2.trialBasis), testBasis_x2(_localoperator2.testBasis),
  J(4),
  hashTableLargeLength(6151), hashTableSmallLength(193)
  //hashTableLargeLength(24593), hashTableSmallLength(769)
{

}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::setJ(int _J)
{
    J = _J;
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::eval(const Coefficients<Lexicographical,T,Index2D> &v,
       Coefficients<Lexicographical,T,Index2D> &AAv)
{
    Coefficients<Lexicographical,T,Index2D> intermediate((size_t)SIZELARGEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index1D> Pe1_UIv(SIZELARGEHASHINDEX1D);

    initializeIntermediateVectorIAv(v, AAv, intermediate);
    evalIA(v, intermediate);
    evalLI(intermediate, AAv);

    intermediate.clear();

    initializeIntermediateVectorUIv(v, AAv, Pe1_UIv);
    evalUI(v, Pe1_UIv, intermediate);
    evalIA(intermediate, AAv);
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::eval(const Coefficients<Lexicographical,T,Index2D> &v,
       Coefficients<Lexicographical,T,Index2D> &AAv,
       T &time_intermediate1, T &time_intermediate2,
       T &time_IAv1, T &time_IAv2, T &time_LIv, T &time_UIv) /*const*/
{
    Coefficients<Lexicographical,T,Index2D> intermediate((size_t)SIZELARGEHASHINDEX2D);

    Coefficients<Lexicographical,T,Index1D> Pe1_UIv(SIZELARGEHASHINDEX1D);

    Timer time;
    time.start();
    initializeIntermediateVectorIAv(v, AAv, intermediate);
    time.stop();
    time_intermediate1 = time.elapsed();

    time.start();
    evalIA(v, intermediate);
    time.stop();
    time_IAv1 = time.elapsed();

    time.start();
    evalLI(intermediate, AAv);
    time.stop();
    time_LIv = time.elapsed();

    intermediate.clear();

    time.start();
    initializeIntermediateVectorUIv(v, AAv, Pe1_UIv);
    time.stop();
    time_intermediate2 = time.elapsed();

    time.start();
    evalUI(v, Pe1_UIv, intermediate);
    time.stop();
    time_UIv = time.elapsed();

    time.start();
    evalIA(intermediate, AAv);
    time.stop();
    time_IAv2 = time.elapsed();
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::debug_eval(const Coefficients<Lexicographical,T,Index2D> &v,
             Coefficients<Lexicographical,T,Index2D> &IAUIv,
             Coefficients<Lexicographical,T,Index2D> &LIIAv,
             const Coefficients<Lexicographical,T,Index2D> &IAv_ref,
             const Coefficients<Lexicographical,T,Index2D> &LIIAv_ref,
             const Coefficients<Lexicographical,T,Index2D> &UIv_ref,
             const Coefficients<Lexicographical,T,Index2D> &IAUIv_ref,
             const Coefficients<Lexicographical,T,Index2D> &AAv_ref) /*const*/
{

    Coefficients<Lexicographical,T,Index2D> intermediate(std::min((size_t)4*IAUIv.size(),
                                                                  (size_t)SIZELARGEHASHINDEX2D));

    Coefficients<Lexicographical,T,Index2D> diff;
    Timer time;

    time.start();
    initializeIntermediateVectorIAv(v, LIIAv, intermediate);
    time.stop();
    std::cerr << "   set up intermediate vector IAv took " << time.elapsed() << " for #v = " << v.size() << std::endl;


    time.start();
    evalIA(v, intermediate);
    time.stop();
    diff = IAv_ref - intermediate;
    std::cerr << "   evalIA took " << time.elapsed() << " for #v = " << v.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();


    time.start();
    evalLI(intermediate, LIIAv);
    time.stop();
    diff = LIIAv_ref - LIIAv;
    std::cerr << "   evalLI took " << time.elapsed() << " for #IAv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    intermediate.clear();

    Coefficients<Lexicographical,T,Index1D> Pe1_UIv(SIZELARGEHASHINDEX1D);
    std::cerr << "   intermediate.size() = " << intermediate.size() << std::endl;
    time.start();
    initializeIntermediateVectorUIv(v, IAUIv, Pe1_UIv);
    time.stop();
    std::cerr << "   set up intermediate vector UIv took " << time.elapsed() << " for #v = " << v.size() << std::endl;


    time.start();
    evalUI(v, Pe1_UIv, intermediate);
    time.stop();
    diff = UIv_ref - intermediate;
    std::cerr << "   evalUI took " << time.elapsed() << " for #v = " << v.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();


    time.start();
    evalIA(intermediate, IAUIv);
    time.stop();
    diff = IAUIv_ref - IAUIv;
    std::cerr << "   evalIA took " << time.elapsed() << " for #UIv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    std::cerr << "   evalIA output size " << IAUIv.size() << std::endl;
    diff.setToZero();

    diff  = AAv_ref - IAUIv;
    diff -= LIIAv;
    std::cerr << "   Diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();
    diff  = AAv_ref - IAUIv_ref;
    diff -= LIIAv_ref;
    std::cerr << "   Diff_ref = " << diff.norm(2.) << std::endl;

}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::initializeIntermediateVectorIAv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &LIIAv,
                                  Coefficients<Lexicographical,T,Index2D> &IAv) const
{
    IAv.clear();
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    XOneAlignedCoefficients x1aligned_LIIAv(n1,n2);
    x1aligned_LIIAv.align(LIIAv,J);




    Coefficients<Lexicographical,T,Index1D> Pe1_v(n1);
    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        if (Pe1_v.find((*col).first.index1)==Pe1_v.end()) {
            Pe1_v[(*col).first.index1] = 0.;
        }
    }

    for (const_coeff1d_it row_x=Pe1_v.begin(); row_x!=Pe1_v.end(); ++row_x) {
        XType xtype_row_x = (*row_x).first.xtype;
        int   j_row_x = (*row_x).first.j;
        long  k_row_x = (*row_x).first.k;
        //IndexSet<Index1D> Lambda_y(98317);
        IndexSet<Index1D> Lambda_y;
        if (xtype_row_x==XWavelet) {
            int j_col_x = 0;
            long k_col_x_first = 0, k_col_x_last = 0;
            trialBasis_x1.getHigherWaveletNeighborsForWavelet(j_row_x, k_row_x, testBasis_x1,
                                                              j_col_x,k_col_x_first,k_col_x_last);
            assert(j_row_x == j_col_x-1);
            Support<T> supp_row_x = trialBasis_x1.psi.support(j_row_x,k_row_x);
            for (long k_col_x=k_col_x_first; k_col_x<=k_col_x_last; ++k_col_x) {
                if (overlap(supp_row_x,testBasis_x1.psi.support(j_col_x,k_col_x))>0) {
                    Index1D col_x(j_col_x,k_col_x,XWavelet);
                    typename XOneAlignedCoefficients::const_map_prindex_it it=x1aligned_LIIAv.map.find(col_x);
                    if (it!=x1aligned_LIIAv.map.end()) {
                        for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                            Lambda_y.insert((*row_y).first);
                        }
                    }
                }
            }
        }
        else {
            int j_col_x = 0;
            long k_col_x_first = 0, k_col_x_last = 0;
            trialBasis_x1.getWaveletNeighborsForScaling(j_row_x, k_row_x, testBasis_x1,
                                                        j_col_x,k_col_x_first,k_col_x_last);
            assert(j_row_x == j_col_x);
            Support<T> supp_row_x = trialBasis_x1.mra.phi.support(j_row_x,k_row_x);
            for (long k_col_x=k_col_x_first; k_col_x<=k_col_x_last; ++k_col_x) {
                if (overlap(supp_row_x,testBasis_x1.psi.support(j_col_x,k_col_x))>0) {
                    Index1D col_x(j_col_x,k_col_x,XWavelet);
                    typename XOneAlignedCoefficients::const_map_prindex_it it=x1aligned_LIIAv.map.find(col_x);
                    if (it!=x1aligned_LIIAv.map.end()) {
                        for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                            Lambda_y.insert((*row_y).first);
                        }
                    }
                }
            }
        }
        for (const_set1d_it row_y=Lambda_y.begin(); row_y!=Lambda_y.end(); ++row_y) {
            IAv[Index2D((*row_x).first,(*row_y))] = 0.;
        }
    }

}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::initializeIntermediateVectorUIv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &IAUIv,
                                  Coefficients<Lexicographical,T,Index1D> &Pe1_UIv) const
{
    for (const_coeff2d_it it=IAUIv.begin(); it!=IAUIv.end(); ++it) {
        Pe1_UIv[(*it).first.index1] = 0.;
    }
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::evalIA(const Coefficients<Lexicographical,T,Index2D> &z,
         Coefficients<Lexicographical,T,Index2D> &IAz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    XOneAlignedCoefficients x1aligned_z(n1,n2);
    XOneAlignedCoefficients x1aligned_IAz(n1,n2);
    x1aligned_z.align(z,J);
    x1aligned_IAz.align(IAz,J);
    time.stop();
    T time_x1align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;

    for (typename XOneAlignedCoefficients::const_map_prindex_it it=x1aligned_z.map.begin();
                                                            it!=x1aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_x = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x2(n2,trialBasis_x2.j0);
        PsiLambdaHat_x2 = (*it).second;
        TreeCoefficients1D<T> PsiLambdaCheck_x2(n2,testBasis_x2.j0);

        if ( x1aligned_IAz.map.find((*it).first)==x1aligned_IAz.map.end() ) continue;
        PsiLambdaCheck_x2 = x1aligned_IAz.map[(*it).first];
        int maxTreeLevel = PsiLambdaCheck_x2.getMaxTreeLevel();
        PsiLambdaCheck_x2.setToZero();
        time.stop();
        time_setup_tree += time.elapsed();


        time.start();
        localoperator2.eval(PsiLambdaHat_x2, PsiLambdaCheck_x2, "A");
        time.stop();
        time_mv1d += time.elapsed();


        time.start();
        PsiLambdaCheck_x2.template addTo<Index2D,Index1D,XOne>(row_x,IAz);
        time.stop();
        time_add_aligned += time.elapsed();
    }

    /*
    std::cerr << "      evalIA: x1align of v took       " << time_x1align_v << std::endl;
    std::cerr << "      evalIA: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalIA: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalIA: add aligned result      " << time_add_aligned << IAz.size() << std::endl;
    */
    return;
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::evalLI(const Coefficients<Lexicographical,T,Index2D> &z,
         Coefficients<Lexicographical,T,Index2D> &LIz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    XTwoAlignedCoefficients x2aligned_z(n1,n2);
    XTwoAlignedCoefficients x2aligned_LIz(n1,n2);
    x2aligned_z.align(z,J);
    x2aligned_LIz.align(LIz,J);
    time.stop();
    T time_x2align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;

    for (typename XTwoAlignedCoefficients::const_map_prindex_it it=x2aligned_z.map.begin();
                                                            it!=x2aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x1(n2,trialBasis_x1.j0);
        PsiLambdaHat_x1 = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();


        time.start();
        TreeCoefficients1D<T> PsiLambdaCheck_x1(n2,testBasis_x1.j0);
        PsiLambdaCheck_x1 = x2aligned_LIz.map[(*it).first];
        int maxTreeLevel = PsiLambdaCheck_x1.getMaxTreeLevel();
        PsiLambdaCheck_x1.setToZero();
        time.stop();

        time.start();
        localoperator1.eval(PsiLambdaHat_x1, PsiLambdaCheck_x1, "L");
        time.stop();
        time_mv1d += time.elapsed();

        time.start();
        PsiLambdaCheck_x1.template addTo<Index2D,Index1D,XTwo>(row_y,LIz);
        time.stop();
        time_add_aligned += time.elapsed();


    }
/*
    std::cerr << "      evalLI: x2align of v took       " << time_x2align_v << std::endl;
    std::cerr << "      evalLI: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalLI: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalLI: add aligned result      " << time_add_aligned << std::endl;
*/
    return;
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<LocalOperator1, LocalOperator2>
::evalUI(const Coefficients<Lexicographical,T,Index2D> &z,
         const Coefficients<Lexicographical,T,Index1D> &Pe1_UIz,
         Coefficients<Lexicographical,T,Index2D> &UIz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    XTwoAlignedCoefficients x2aligned_z(n1,n2);
    x2aligned_z.align(z,J);
    time.stop();
    T time_x2align_z = time.elapsed();

    T time_initial_outputset = 0.;
    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;


    for (typename XTwoAlignedCoefficients::const_map_prindex_it it=x2aligned_z.map.begin();
                                                            it!=x2aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x1(n2,trialBasis_x1.j0);
        PsiLambdaHat_x1 = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        int maxTreeLevel = PsiLambdaHat_x1.getMaxTreeLevel();


        time.start();
        TreeCoefficients1D<T> PsiLambdaCheck_x1(n2,testBasis_x1.j0);
        // Checking scaling functions
        for (const_by_level_it level_it =PsiLambdaHat_x1[0].map.begin();
                               level_it!=PsiLambdaHat_x1[0].map.end(); ++level_it) {
            int  j_scaling1 = trialBasis_x1.j0;
            long k_scaling1 = (*level_it).first;
            int  j_scaling2 = 0;
            long k_scaling_first = 0, k_scaling_last = 0;
            trialBasis_x1.getScalingNeighborsForScaling(j_scaling1,k_scaling1, testBasis_x1,
                                                        j_scaling2,k_scaling_first,k_scaling_last);
            assert(j_scaling1==j_scaling2);
            Support<T> supp1 = trialBasis_x1.mra.phi.support(j_scaling1,k_scaling1);
            for (int k_scaling2=k_scaling_first; k_scaling2<=k_scaling_last; ++k_scaling2) {
                if (   overlap(supp1, testBasis_x1.mra.phi.support(j_scaling2,k_scaling2) ) >0
                    && Pe1_UIz.find(Index1D(j_scaling2,k_scaling2,XBSpline))!=Pe1_UIz.end() ) {
                    PsiLambdaCheck_x1[0].map[k_scaling2] = 0.;
                }
            }
        }
        for (int i=1; i<=maxTreeLevel; ++i) {
            for (const_by_level_it level_it =PsiLambdaHat_x1[i].map.begin();
                                   level_it!=PsiLambdaHat_x1[i].map.end(); ++level_it) {
                int  j_wavelet1 = trialBasis_x1.j0+i-1;
                long k_wavelet1 = (*level_it).first;
                int  j_wavelet2 = 0;
                long k_wavelet_first = 0, k_wavelet_last = 0;
                trialBasis_x1.getWaveletNeighborsForWavelet(j_wavelet1,k_wavelet1, testBasis_x1,
                                                            j_wavelet2,k_wavelet_first,k_wavelet_last);
                assert(j_wavelet1==j_wavelet2);
                Support<T> supp1 = trialBasis_x1.psi.support(j_wavelet1,k_wavelet1);
                for (int k_wavelet2=k_wavelet_first; k_wavelet2<=k_wavelet_last; ++k_wavelet2) {
                    if (   overlap(supp1,testBasis_x1.psi.support(j_wavelet2,k_wavelet2) ) >0
                        && Pe1_UIz.find(Index1D(j_wavelet2,k_wavelet2,XWavelet))!=Pe1_UIz.end() ) {
                        PsiLambdaCheck_x1[i].map[k_wavelet2] = 0.;
                    }
                }
            }
        }

        PsiLambdaCheck_x1.setMaxTreeLevel(maxTreeLevel);
        time.stop();
        time_initial_outputset += time.elapsed();

        time.start();
        localoperator1.eval(PsiLambdaHat_x1, PsiLambdaCheck_x1, "U");
        time.stop();
        time_mv1d += time.elapsed();


        time.start();
        PsiLambdaCheck_x1.template addTo<Index2D,Index1D,XTwo>(row_y,UIz);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    /*
    std::cerr << "      evalUI: x2align of v took       " << time_x2align_z << std::endl;
    std::cerr << "      evalUI: time for initial output " << time_initial_outputset << std::endl;
    std::cerr << "      evalUI: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalUI: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalUI: add aligned result      " << time_add_aligned << std::endl;
    */
    return;
}

}   // namespace lawa


/*
    int count = 0;
    time.start();
    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_LIIAv.map.begin(); it!=x1aligned_LIIAv.map.end(); ++it) {
        Index1D col_x = (*it).first;
        if (col_x.xtype==XWavelet && col_x.j>trialBasis_x1.j0) {
            int j_row = 0;
            long k_row_first = 0, k_row_last = 0;
            testBasis_x1.getLowerWaveletNeighborsForWavelet(col_x.j, col_x.k, trialBasis_x1,j_row,k_row_first,k_row_last);
            assert(j_row == col_x.j-1);
            Support<T> supp_col_x = testBasis_x1.psi.support(col_x.j,col_x.k);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {

                Index1D row_x(j_row,k_row,XWavelet);
                if (  (Pe1_v.find(row_x)!=Pe1_v.end()) &&  overlap(trialBasis_x1.psi.support(row_x.j,row_x.k), supp_col_x)>0
                    //&& (Tree_Pe1_v[j_row-trialBasis_x1.j0+1].map.find(k_row)!=Tree_Pe1_v[j_row-trialBasis_x1.j0+1].map.end() ) ) {
                    ) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        Index2D index(row_x,(*row_y).first);
                        if (IAv.find(index) == IAv.end()) {
                            IAv[Index2D(row_x,(*row_y).first)] = 0.;
                        }
                        ++count;
                    }
                }
            }
        }
        else if (col_x.xtype==XWavelet && col_x.j==trialBasis_x1.j0) {
            int j_row = 0;
            long k_row_first = 0, k_row_last = 0;
            testBasis_x1.getScalingNeighborsForWavelet(col_x.j,col_x.k,trialBasis_x1,j_row,
                                                                    k_row_first,k_row_last);
            assert(j_row == col_x.j);
            Support<T> supp_col_x = testBasis_x1.psi.support(col_x.j,col_x.k);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                Index1D row_x(j_row,k_row,XBSpline);
                if ( Pe1_v.find(row_x)!=Pe1_v.end() && overlap(trialBasis_x1.mra.phi.support(row_x.j,row_x.k), supp_col_x)>0
                    //&& (Tree_Pe1_v[0].map.find(k_row)!=Tree_Pe1_v[0].map.end() ) ) {
                    ) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        Index2D index(row_x,(*row_y).first);
                        if (IAv.find(index) == IAv.end()) {
                            IAv[Index2D(row_x,(*row_y).first)] = 0.;
                        }

                        ++count;
                    }
                }
            }
        }
    }
    time.stop();
    std::cerr << "   -> Set up: " << time.elapsed() << ", #IAv = " << IAv.size() << " " << count << std::endl;
    */
/*
    time.start();
    alignedCoefficients2 x1aligned_LIIAv2;
    x1aligned_LIIAv2.align_x1(LIIAv);
    time.stop();
    std::cerr << "   -> New Alignment: " << time.elapsed() << std::endl;
    time.start();
    int count = 0;
    for (typename alignedCoefficients2::const_coeff_prinindex_it it=x1aligned_LIIAv2.principalIndices.begin();
                 it!=x1aligned_LIIAv2.principalIndices.end(); ++it) {
       Index1D col_x = (*it).first;
       int pos = (*it).second;
       if (col_x.xtype==XWavelet && col_x.j>trialBasis_x1.j0) {
           int j_row = 0;
           long k_row_first = 0, k_row_last = 0;
           testBasis_x1.getLowerWaveletNeighborsForWavelet(col_x.j, col_x.k, trialBasis_x1,j_row,k_row_first,k_row_last);
           assert(j_row == col_x.j-1);
           Support<T> supp_col_x = testBasis_x1.psi.support(col_x.j,col_x.k);
           for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {

               Index1D row_x(j_row,k_row,XWavelet);
               if (  (Pe1_v.find(row_x)!=Pe1_v.end()) &&  overlap(trialBasis_x1.psi.support(row_x.j,row_x.k), supp_col_x)>0
                   //&& (Tree_Pe1_v[j_row-trialBasis_x1.j0+1].map.find(k_row)!=Tree_Pe1_v[j_row-trialBasis_x1.j0+1].map.end() ) ) {
                   ) {
                   //typename alignedCoefficients2::AlignedIndices *p_alignedIndices = &x1aligned_LIIAv2.principalIndexToAlignedIndices[pos];
                   for (typename alignedCoefficients2::AlignedIndices::const_iterator  aligned_it=x1aligned_LIIAv2.principalIndexToAlignedIndices[pos].begin();
                           aligned_it != x1aligned_LIIAv2.principalIndexToAlignedIndices[pos].end(); ++aligned_it) {
                       Index1D row_y = *(*aligned_it);
                       //IAv[Index2D(row_x,row_y)] = 0.;
                       ++count;
                   }
               }
           }
       }
       else if (col_x.xtype==XWavelet && col_x.j==trialBasis_x1.j0) {
           int j_row = 0;
           long k_row_first = 0, k_row_last = 0;
           testBasis_x1.getScalingNeighborsForWavelet(col_x.j,col_x.k,trialBasis_x1,j_row,
                                                                   k_row_first,k_row_last);
           assert(j_row == col_x.j);
           Support<T> supp_col_x = testBasis_x1.psi.support(col_x.j,col_x.k);
           for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
               Index1D row_x(j_row,k_row,XBSpline);
               if ( Pe1_v.find(row_x)!=Pe1_v.end() && overlap(trialBasis_x1.mra.phi.support(row_x.j,row_x.k), supp_col_x)>0
                   //&& (Tree_Pe1_v[0].map.find(k_row)!=Tree_Pe1_v[0].map.end() ) ) {
                   ) {
                   //typename alignedCoefficients2::AlignedIndices *p_alignedIndices = &x1aligned_LIIAv2.principalIndexToAlignedIndices[pos];
                   for (typename alignedCoefficients2::AlignedIndices::const_iterator  aligned_it=x1aligned_LIIAv2.principalIndexToAlignedIndices[pos].begin();
                          aligned_it != x1aligned_LIIAv2.principalIndexToAlignedIndices[pos].end(); ++aligned_it) {
                       Index1D row_y = *(*aligned_it);
                       //IAv[Index2D(row_x,row_y)] = 0.;
                       ++count;
                   }
               }
           }
       }
   }
   time.stop();
   std::cerr << "   -> Set up: " << time.elapsed() << " " << count << std::endl;

*/
