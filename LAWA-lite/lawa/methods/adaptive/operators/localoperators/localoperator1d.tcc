namespace lawa {

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm,
          typename BilinearForm>
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::LocalOperator1D(const TestBasis &_testBasis, const TrialBasis &_trialBasis,
                  RefinementBilinearForm &_RefinementBil)
: testBasis(_testBasis), trialBasis(_trialBasis),
  testRefinementBasis(testBasis.refinementbasis), trialRefinementBasis(trialBasis.refinementbasis),
  RefinementBil(_RefinementBil), Bil(_RefinementBil),
  testLocalRefine(testBasis), trialLocalRefine(trialBasis),
  testRefinementLevelOffset(testBasis.j0), trialRefinementLevelOffset(trialBasis.j0)
{
    assert(testBasis.j0==trialBasis.j0);    // Using different scales is not meaningful.
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm,
          typename BilinearForm>
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::LocalOperator1D(const TestBasis &_testBasis, const TrialBasis &_trialBasis,
                  RefinementBilinearForm &_RefinementBil, BilinearForm &_Bil)
: testBasis(_testBasis), trialBasis(_trialBasis),
  testRefinementBasis(testBasis.refinementbasis), trialRefinementBasis(trialBasis.refinementbasis),
  RefinementBil(_RefinementBil), Bil(_Bil),
  testLocalRefine(testBasis), trialLocalRefine(trialBasis),
  testRefinementLevelOffset(testBasis.j0), trialRefinementLevelOffset(trialBasis.j0)
{
    assert(testBasis.j0==trialBasis.j0);    // Using different scales is not meaningful.
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm,
          typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::eval(const TreeCoefficients1D<T> &PsiLambdaHat, TreeCoefficients1D<T> &PsiLambdaCheck,
       const char* mode)
{
    CoefficientsByLevel<T> d, PhiPiCheck;
    if (TrialBasis::Cons==Multi) {
        int j_bspline = 0;
        if (PsiLambdaHat[0].map.size()==0) {
            std::cerr << "ERROR: No scaling functions in tree PsiLambdaHat!" << std::endl;
            std::cerr << PsiLambdaHat << std::endl;
        }
        trialLocalRefine.reconstructOnlyMultiScaling(PsiLambdaHat[0], trialBasis.j0, d, j_bspline);
        trialRefinementLevelOffset = j_bspline; //required when using multiwavelets!
    }
    else {
        d = PsiLambdaHat[0];
    }
    if (TestBasis::Cons==Multi) {
        int j_bspline = 0;
        if (strcmp(mode,"L")!=0) {
            if (PsiLambdaCheck[0].map.size()==0) {
                std::cerr << "ERROR: No scaling functions in tree PsiLambdaCheck!" << std::endl;
                std::cerr << PsiLambdaCheck << std::endl;
            }
            testLocalRefine.reconstructOnlyMultiScaling(PsiLambdaCheck[0], testBasis.j0,
                                                PhiPiCheck, j_bspline);
            testRefinementLevelOffset = j_bspline; //required when using multiwavelets!
        }
        else {
            testRefinementLevelOffset = testBasis.mra.phi.getRefinementLevel(testBasis.j0);
        }
    }
    else {
        PhiPiCheck = PsiLambdaCheck[0];
    }

    if (strcmp(mode,"A")==0)                this->_evalA(0, d, PsiLambdaHat, PhiPiCheck, PsiLambdaCheck);
    //else if (strcmp(mode,"A_nonRec")==0)    this->_evalA_nonRecursive(d, PsiLambdaHat, PhiPiCheck, PsiLambdaCheck);
    else if (strcmp(mode,"U")==0)           this->_evalU(0, d, PsiLambdaHat, PhiPiCheck, PsiLambdaCheck);
    else if (strcmp(mode,"L")==0)           this->_evalL(0, d, PsiLambdaHat, PsiLambdaCheck);
    else {
        std::cerr << "LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm>: Unknow mode"
                  << mode << ". Exit." << std::endl;
        exit(1);
        return;
    }

    if (TestBasis::Cons==Multi && strcmp(mode,"L")!=0) {
        testLocalRefine.decompose_OnlyMultiScaling(PhiPiCheck, PsiLambdaCheck[0], testBasis.j0);
    }
    else {
        PsiLambdaCheck[0] = PhiPiCheck;
    }
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_evalA(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck)
{
    Timer time;
    int shifted_l = l+1; // by convention v[0] contains scaling / B-spline indices.
    if (d.map.size()==0 && c[shifted_l].map.size()==0) return;
    if (PhiPiCheck.map.size()==0 && PsiLambdaCheck[shifted_l].map.size()==0) return;

    int j_bspline_test =  l + testRefinementLevelOffset;
    int j_bspline_trial = l + trialRefinementLevelOffset;
    int j_wavelet_test =  l + testBasis.j0;
    int j_wavelet_trial = l + trialBasis.j0;

    //size_t hm_size = COEFFBYLEVELSIZE;
    size_t hm_size = (size_t)(2*std::max(d.map.size(), PhiPiCheck.map.size()));

    // Splitting of B-spline index set $\check{\Phi}$.
    CoefficientsByLevel<T> PhiPiCheck2(l,hm_size);

    /// Optimization: does not compile / work for realline bases
    //if (c[shifted_l].map.size()==trialBasis.cardJ(j_wavelet_trial)) {
    //    PhiPiCheck2.map.swap(PhiPiCheck.map);
    //}
    //else {
    _splitPhiPi(l, c[shifted_l], PhiPiCheck, PhiPiCheck2);
    //}

    // Compute underlinePiCheck
    int j_refinement_test = 0;
    CoefficientsByLevel<T> PhiPiunderlineCheck(l,hm_size);
    testLocalRefine.reconstruct(PhiPiCheck2,         j_bspline_test,
                                PsiLambdaCheck[shifted_l],   j_wavelet_test,
                                PhiPiunderlineCheck, j_refinement_test);

    //Compute a(\check{Phi}|_{\check{Pi}^{(1)}} , \hat{\Phi}|_{\hat{Pi}}) d
    //time.start();
    _applyRefinementBilinearForm(l, d, PhiPiCheck);
    //time.stop();

    // Splitting of B-spline coefficient vector $d$.
    //time.start();
    CoefficientsByLevel<T> d2(l,hm_size);
    /// Optimization: does not compile / work for realline bases
    //if (PsiLambdaCheck[shifted_l].map.size()==testBasis.cardJ(j_wavelet_test)) {
    //    d2.map.swap(d.map);
    //}
    //else {
    _splitd(l, PsiLambdaCheck[shifted_l], d, d2);
    //}
    //time.stop();
    //std::cerr << "l = " << l << " : splitting of d took " << time.elapsed() << std::endl;

    // Compute underline d
    //time.start();
    CoefficientsByLevel<T> underline_d(l+1,4*hm_size);
    int j_refinement_trial = 0;
    trialLocalRefine.reconstruct(d2,           j_bspline_trial,
                                 c[shifted_l], j_wavelet_trial,
                                 underline_d,  j_refinement_trial);
    //time.stop();
    //std::cerr << "l = " << l << " : trialLocalRefine.reconstruct took " << time.elapsed()
    //          << ", output size: " << underline_d.map.size() << std::endl;

    //Recursive call of eval
    this->_evalA(l+1, underline_d, c, PhiPiunderlineCheck, PsiLambdaCheck);

    testLocalRefine.decompose_(PhiPiunderlineCheck,
                               PhiPiCheck2,               j_bspline_test,
                               PsiLambdaCheck[shifted_l], j_wavelet_test);

    //Compute a(\check{Phi}|_{\check{Pi}^{(2)}} , \hat{\Phi}|_{\hat{Pi}^{(1)}}) d^{(1)}
    _applyRefinementBilinearForm(l, d, PhiPiCheck2);

    /// Optimization: does not compile / work for realline bases
    //if (c[shifted_l].map.size()==trialBasis.cardJ(j_wavelet_trial)) {
    //    PhiPiCheck.map.swap(PhiPiCheck2.map);
    //}
    //else {
    PhiPiCheck += PhiPiCheck2;
    //}

    return;

}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_evalU(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck)
{
    int shifted_l = l+1; // by convention v[0] contains scaling / B-spline indices.
    if (d.map.size()==0 && c[shifted_l].map.size()==0) return;

    int j_bspline_test =  l + testRefinementLevelOffset;
    int j_bspline_trial = l + trialRefinementLevelOffset;
    int j_wavelet_test =  l + testBasis.j0;
    int j_wavelet_trial = l + trialBasis.j0;

    size_t hm_size = 2*COEFFBYLEVELSIZE;

    // Splitting of B-spline index set $\check{\Phi}$.
    CoefficientsByLevel<T> PhiPiCheck2(l,hm_size);
    _splitPhiPi(l, c[shifted_l], PhiPiCheck, PhiPiCheck2);

    // Compute underlinePiCheck
    int j_refinement_test = 0;
    CoefficientsByLevel<T> PhiPiunderlineCheck(l,hm_size);
    testLocalRefine.reconstruct(PhiPiCheck2,         j_bspline_test,
                                PsiLambdaCheck[shifted_l],   j_wavelet_test,
                                PhiPiunderlineCheck, j_refinement_test);

    //Compute a(\check{Phi}|_{\check{Pi}^{(1)}} , \hat{\Phi}|_{\hat{Pi}}) d
    _applyRefinementBilinearForm(l, d, PhiPiCheck);

    // Compute underline d
    CoefficientsByLevel<T> help, underline_d(l+1,hm_size);
    int j_refinement_trial = 0;
    trialLocalRefine.reconstruct(help,         j_bspline_trial,
                                 c[shifted_l], j_wavelet_trial,
                                 underline_d,  j_refinement_trial);

    //Recursive call of eval
    this->_evalU(l+1, underline_d, c, PhiPiunderlineCheck, PsiLambdaCheck);

    testLocalRefine.decompose_(PhiPiunderlineCheck,
                               PhiPiCheck2,               j_bspline_test,
                               PsiLambdaCheck[shifted_l], j_wavelet_test);

    //Compute a(\check{Phi}|_{\check{Pi}^{(2)}} , \hat{\Phi}|_{\hat{Pi}^{(1)}}) d^{(1)}
    _applyRefinementBilinearForm(l, d, PhiPiCheck2);

    PhiPiCheck += PhiPiCheck2;

    return;
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_evalL(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
         TreeCoefficients1D<T> &PsiLambdaCheck)
{
    int shifted_l = l+1; // by convention v[0] contains scaling / B-spline indices.
    if (d.map.size()==0 && c[shifted_l].map.size()==0) return;

    int j_bspline_test =  l + testRefinementLevelOffset;
    int j_bspline_trial = l + trialRefinementLevelOffset;
    int j_wavelet_test =  l + testBasis.j0;
    int j_wavelet_trial = l + trialBasis.j0;

    size_t hm_size = 2*COEFFBYLEVELSIZE;

    // Compute underlinePiCheck
    int j_refinement_test = 0;
    CoefficientsByLevel<T> help, PhiPiunderlineCheck(l,hm_size);
    testLocalRefine.reconstruct(help,                      j_bspline_test,
                                PsiLambdaCheck[shifted_l], j_wavelet_test,
                                PhiPiunderlineCheck,       j_refinement_test);

    // Splitting of B-spline coefficient vector $d$.
    CoefficientsByLevel<T> d2(l,hm_size);
    _splitd(l, PsiLambdaCheck[shifted_l], d, d2);

    // Compute $\underline{\underline d}$
    int j_refinement_trial = 0;
    CoefficientsByLevel<T> underline_d;
    trialLocalRefine.reconstruct(d2,          j_bspline_trial,
                                 help,        j_wavelet_trial,
                                 underline_d, j_refinement_trial);

    _applyRefinementBilinearForm(l+1, underline_d, PhiPiunderlineCheck);

    testLocalRefine.decompose_(PhiPiunderlineCheck,
                               help,                      j_bspline_test,
                               PsiLambdaCheck[shifted_l], j_wavelet_test);

    // Compute underline d
    trialLocalRefine.reconstruct(help,         j_bspline_trial,
                                 c[shifted_l], j_wavelet_trial,
                                 underline_d,  j_refinement_trial);

    //Recursive call of eval
    this->_evalL(l+1, underline_d, c, PsiLambdaCheck);

    return;
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_splitPhiPi(int l, const CoefficientsByLevel<T> &c_l, CoefficientsByLevel<T> &PhiPiCheck1,
                 CoefficientsByLevel<T> &PhiPiCheck2) const
{
    if (c_l.map.size()==0) return;

    if (c_l.map.size()>PhiPiCheck1.map.size()) {
        const_by_level_it p_c_l_end = c_l.map.end();
        int j_bspline = l + testRefinementLevelOffset;
        for (const_by_level_it mu=PhiPiCheck1.map.begin(); mu!=PhiPiCheck1.map.end(); ++mu) {
            long k_bspline = (*mu).first;
            int  j_wavelet = 0;
            long k_wavelet_first = 0L, k_wavelet_last = 0L;
            testRefinementBasis.getWaveletNeighborsForBSpline(j_bspline, k_bspline, trialBasis,
                                                    j_wavelet, k_wavelet_first, k_wavelet_last);
            for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
                const_by_level_it p_entry_c_l = c_l.map.find(k_wavelet);
                if (p_entry_c_l!=p_c_l_end) {
                    PhiPiCheck2.map[k_bspline] = 0.;
                    PhiPiCheck1.map.erase(k_bspline);
                    break;
                }
            }
        }
    }
    else {
        int j_wavelet = l + trialBasis.j0;
        for (const_by_level_it mu=c_l.map.begin(); mu!=c_l.map.end(); ++mu) {
            long k_wavelet = (*mu).first;
            int  j_bspline = 0;
            long k_bspline_first = 0L, k_bspline_last = 0L;
            trialBasis.getBSplineNeighborsForWavelet(j_wavelet, k_wavelet, testRefinementBasis,
                                                     j_bspline, k_bspline_first, k_bspline_last);
            for (long k_bspline=k_bspline_first; k_bspline<=k_bspline_last; ++k_bspline) {
                const_by_level_it p_PhiPiCheck1 = PhiPiCheck1.map.find(k_bspline);
                if (p_PhiPiCheck1!=PhiPiCheck1.map.end()) {
                    PhiPiCheck2.map[k_bspline] = 0.;
                    PhiPiCheck1.map.erase(k_bspline);
                }
            }
        }
    }
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_splitd(int l, const CoefficientsByLevel<T> &PsiLambdaCheck_l,
          CoefficientsByLevel<T> &d1, CoefficientsByLevel<T> &d2) const
{
    if (PsiLambdaCheck_l.map.size()==0) return;
    if (d1.map.size()>PsiLambdaCheck_l.map.size()) {
        int j_wavelet = l + testBasis.j0;
        for (const_by_level_it mu=PsiLambdaCheck_l.map.begin(); mu!=PsiLambdaCheck_l.map.end(); ++mu) {
            long k_wavelet = (*mu).first;
            int  j_bspline = 0;
            long k_bspline_first = 0L, k_bspline_last = 0L;
            testBasis.getBSplineNeighborsForWavelet(j_wavelet, k_wavelet, trialRefinementBasis,
                                                    j_bspline, k_bspline_first, k_bspline_last);
            //Support<T> testSupp = testBasis.generator(XWavelet).support(j_wavelet,k_wavelet);
            for (long k_bspline=k_bspline_first; k_bspline<=k_bspline_last; ++k_bspline) {
                //Support<T> trialSupp = trialBasis.generator(XBSpline).support(j_bspline,k_bspline);
                //if (!overlap(testSupp,trialSupp)) continue;
                by_level_it p_entry_d1 = d1.map.find(k_bspline);
                if (p_entry_d1!=d1.map.end()) {
                    d2.map[k_bspline] = (*p_entry_d1).second;
                    d1.map.erase(p_entry_d1);
                }
            }
        }
    }
    else {
        const_by_level_it p_PsiLambdaCheck_l_end = PsiLambdaCheck_l.map.end();
        int j_bspline = l + trialRefinementLevelOffset;
        for (const_by_level_it mu=d1.map.begin(); mu!=d1.map.end(); ++mu) {
            long k_bspline = (*mu).first;
            int  j_wavelet = 0;
            long k_wavelet_first = 0L, k_wavelet_last = 0L;
            trialRefinementBasis.getWaveletNeighborsForBSpline(j_bspline, k_bspline, testBasis,
                                                     j_wavelet, k_wavelet_first, k_wavelet_last);
            //Support<T> trialSupp = trialBasis.generator(XBSpline).support(j_bspline,k_bspline);
            for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
                //Support<T> testSupp = testBasis.generator(XWavelet).support(j_wavelet,k_wavelet);
                const_by_level_it p_PsiLambdaCheck_l = PsiLambdaCheck_l.map.find(k_wavelet);
                if (/*overlap(testSupp,trialSupp) &&*/ p_PsiLambdaCheck_l!=p_PsiLambdaCheck_l_end) {
                    d2.map[k_bspline] = (*mu).second;
                    d1.map.erase(k_bspline);
                    break;
                }
            }
        }
    }
}

template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_applyRefinementBilinearForm(int l, const CoefficientsByLevel<T> &d, CoefficientsByLevel<T> &PhiPiCheck)
{
    if (d.map.size()==0 || PhiPiCheck.map.size()==0) return;
    if (d.map.size() > PhiPiCheck.map.size()) {
        const_by_level_it p_d_end = d.map.end();
        int j_bspline1 = l + testRefinementLevelOffset;
        int j_bspline2=0;
        long k_bspline2_first=0L, k_bspline2_last=0L;

        for (by_level_it mu=PhiPiCheck.map.begin(); mu!=PhiPiCheck.map.end(); ++mu) {
            long double val = 0.;
            long k_bspline1 = (*mu).first;

            testRefinementBasis.getBSplineNeighborsForBSpline(j_bspline1, k_bspline1,
                               trialRefinementBasis, j_bspline2, k_bspline2_first, k_bspline2_last);

            for (long k_bspline2=k_bspline2_first; k_bspline2<=k_bspline2_last; ++k_bspline2) {
                const_by_level_it p_entry_d = d.map.find(k_bspline2);
                if (p_entry_d!=p_d_end) {
                    val += (long double)(RefinementBil(XBSpline, j_bspline1, k_bspline1,
                                                       XBSpline, j_bspline2, k_bspline2) * (*p_entry_d).second);
                }
            }
            (*mu).second += val;
        }
    }

    else {
        //const_by_level_it p_PhiPiCheck_end = PhiPiCheck.map.end();
        int j_bspline2 = l + trialRefinementLevelOffset;
        int j_bspline1=0;
        long k_bspline1_first=0L, k_bspline1_last=0L;

        for (const_by_level_it mu=d.map.begin(); mu!=d.map.end(); ++mu) {
            long double val = 0.;
            long k_bspline2 = (*mu).first;
            trialRefinementBasis.getBSplineNeighborsForBSpline(j_bspline2, k_bspline2,
                                testRefinementBasis, j_bspline1, k_bspline1_first, k_bspline1_last);
            for (long k_bspline1=k_bspline1_first; k_bspline1<=k_bspline1_last; ++k_bspline1) {
                by_level_it p_PhiPiCheck=PhiPiCheck.map.find(k_bspline1);
                if (p_PhiPiCheck!=PhiPiCheck.map.end()) {
                    (*p_PhiPiCheck).second
                    += (long double)(RefinementBil(XBSpline, j_bspline1, k_bspline1,
                                                   XBSpline, j_bspline2, k_bspline2) * (*mu).second);
                }
            }
        }
    }
    /*
    for (by_level_it mu=PhiPiCheck.map.begin(); mu!=PhiPiCheck.map.end(); ++mu) {
        for (const_by_level_it lambda=d.map.begin(); lambda!=d.map.end(); ++lambda) {
            (*mu).second += (long double)(Bil(XBSpline, l + testRefinementLevelOffset, (*mu).first,
                                              XBSpline, l + trialRefinementLevelOffset,(*lambda).first)
                            * (*lambda).second);
        }
    }
    */
}


/*
template <typename TestBasis, typename TrialBasis, typename RefinementBilinearForm, typename BilinearForm>
void
LocalOperator1D<TestBasis, TrialBasis, RefinementBilinearForm, BilinearForm>
::_evalA_nonRecursive(CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
                      CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck)
{
    int lmax = 0;
    TreeCoefficients1D<T> PhiPiCheck1(COEFFBYLEVELSIZE);
    TreeCoefficients1D<T> PhiPiCheck2(COEFFBYLEVELSIZE);
    for (int l=0; l<=100; ++l) {
        //std::cout << "******* l = " << l <<  " *********" << std::endl;
        int shifted_l = l+1; // by convention v[0] contains scaling / B-spline indices.

        if (d.map.size()==0 && c[shifted_l].map.size()==0) break;

        int j_bspline_test =  l + testRefinementLevelOffset;
        int j_bspline_trial = l + trialRefinementLevelOffset;
        int j_wavelet_test =  l + testBasis.j0;
        int j_wavelet_trial = l + trialBasis.j0;

        size_t hm_size = l > 7 ? COEFFBYLEVELSIZE : 255;

        // Splitting of B-spline index set $\check{\Phi}$.
        _splitPhiPi(l, c[shifted_l], PhiPiCheck, PhiPiCheck2[l]);
        PhiPiCheck1[l] = PhiPiCheck;
        //std::cout << " PhiPiCheck1[l] = " << PhiPiCheck1[l] << std::endl;
        //std::cout << " PhiPiCheck2[l] = " << PhiPiCheck2[l] << std::endl;

        // Compute underlinePiCheck
        int j_refinement_test = 0;
        PhiPiCheck.map.clear();
        testLocalRefine.reconstruct(PhiPiCheck2[l],             j_bspline_test,
                                    PsiLambdaCheck[shifted_l],  j_wavelet_test,
                                    PhiPiCheck,                 j_refinement_test);

        //Compute a(\check{Phi}|_{\check{Pi}^{(1)}} , \hat{\Phi}|_{\hat{Pi}}) d
        _applyRefinementBilinearForm(l, d, PhiPiCheck1[l]);

        // Splitting of B-spline coefficient vector $d$.
        CoefficientsByLevel<T> d2(l,hm_size);
        _splitd(l, PsiLambdaCheck[shifted_l], d, d2);
        //std::cout << " d1  = " << d << std::endl;
        //std::cout << " d2  = " << d2 << std::endl;

        //Compute a(\check{Phi}|_{\check{Pi}^{(2)}} , \hat{\Phi}|_{\hat{Pi}^{(1)}}) d^{(1)}
        _applyRefinementBilinearForm(l, d, PhiPiCheck2[l]);

        // Compute underline d
        d.map.clear();
        d.set(l+1,hm_size);
        int j_refinement_trial = 0;
        trialLocalRefine.reconstruct(d2,           j_bspline_trial,
                                     c[shifted_l], j_wavelet_trial,
                                     d,            j_refinement_trial);

        ++lmax;
    }
    //std::cout << std::endl << std::endl;
    //std::cout << "  Going back again..." << std::endl;
    for (int l=lmax+1; l>=1; --l) {
        std::cout << "******* l = " << l <<  " *********" << std::endl;
        int j_bspline_test =  l + testRefinementLevelOffset - 1;
        int j_wavelet_test =  l + testBasis.j0 - 1;

        //std::cout << "PhiPiCheck1[l] = " << PhiPiCheck1[l] << std::endl;
        testLocalRefine.decompose_(PhiPiCheck1[l],
                                   PhiPiCheck2[l-1],            j_bspline_test,
                                   PsiLambdaCheck[l], j_wavelet_test);
        //std::cout << "PhiPiCheck2[l-1]         = " << PhiPiCheck2[l-1] << std::endl;
        //std::cout << "PsiLambdaCheck[l-1]    = " << PsiLambdaCheck[l] << std::endl;

        PhiPiCheck1[l-1] += PhiPiCheck2[l-1];
    }
    PhiPiCheck = PhiPiCheck1[0];
}
*/

}   // namespace lawa
