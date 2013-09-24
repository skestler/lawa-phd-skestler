namespace lawa {

template <typename PrimalBasis>
LocalRefinement<PrimalBasis>::LocalRefinement(const PrimalBasis &_basis)
 : basis(_basis), refinementbasis(basis.refinementbasis)
{

}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T,Index1D> &u, int j_scaling,
                                           Coefficients<Lexicographical,T,Index1D> &u_loc_single) const
{
    TreeCoefficients1D<T> u_tree(255,basis.j0);
    fromCoefficientsToTreeCoefficients(u, u_tree);
    int j_bspline = j_scaling;
    int j_wavelet = j_scaling;

    CoefficientsByLevel<T> u_bspline;
    if (PrimalBasis::Cons==Multi && basis.d>1) {
        this->reconstructOnlyMultiScaling(u_tree.bylevel[0], j_scaling, u_bspline, j_bspline);
        u_tree.bylevel[0] = u_bspline;
    }

    int imax = u_tree.getMaxTreeLevel();
    for (int i=0; i<imax; ++i) {
        int  j_refinement = j_bspline + i;
	CoefficientsByLevel<T> help;
	help = u_tree[i];
        for (const_coeffbylevel_it it=help.map.begin(); it!=help.map.end(); ++it) {
            long k_refinement = (*it).first;
            int test_j_wavelet = 0;
            long k_first = 0L, k_last = 0L;
            refinementbasis.getWaveletNeighborsForBSpline(j_refinement,k_refinement, basis, test_j_wavelet, k_first, k_last);
            assert(test_j_wavelet==j_wavelet+i);
            bool has_neighbor=false;
            for (long k=k_first; k<=k_last; ++k) {
                if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                    has_neighbor = true;
                    break;
                }
            }
            if (!has_neighbor) {
                u_tree[i].map.erase((*it).first);
                u_loc_single[Index1D(j_refinement,k_refinement,XBSpline)] = (*it).second;
            }
        }

        CoefficientsByLevel<T> u_loc_single_jP1;
        int test_j_refinement = 0;
        this->reconstruct(u_tree[i], j_bspline+i, u_tree[i+1], j_wavelet+i, u_loc_single_jP1, test_j_refinement);
        assert(test_j_refinement==j_refinement+1);
        u_tree[i+1] = u_loc_single_jP1;
    }
    for (const_coeffbylevel_it it=u_tree[imax].map.begin(); it!=u_tree[imax].map.end(); ++it) {
        u_loc_single[Index1D(j_bspline+imax,(*it).first,XBSpline)] = (*it).second;
    }
    return;
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const CoefficientsByLevel<T> &u_bspline, int j_bspline,
                                           const CoefficientsByLevel<T> &u_wavelet, int j_wavelet,
                                           CoefficientsByLevel<T> &u_loc_single, int &j_refinement) const
{
    int j1_refinement = refinementbasis.mra.phi.getRefinementLevel(j_bspline);
    // pre-initialization need if we do not enter the following loop.
    for (typename CoefficientsByLevel<T>::const_it it=u_bspline.map.begin(); it!=u_bspline.map.end(); ++it) {
        this->reconstructBSpline(j_bspline, (*it).first, (*it).second, u_loc_single, j1_refinement);
    }
    int j2_refinement = basis.psi.getRefinementLevel(j_wavelet);
    // pre-initialization need if we do not enter the following loop.
    for (typename CoefficientsByLevel<T>::const_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        this->reconstructWavelet(j_wavelet, (*it).first, (*it).second, u_loc_single, j2_refinement);
    }
    //assert(j1_refinement==j2_refinement);
    if(j1_refinement!=j2_refinement) {
        std::cerr << "LocalRefinement<PrimalBasis> ERROR: j_bspline = " << j_bspline << ", j_wavelet = " << j_wavelet << std::endl;
        std::cerr << "                                j1_refinement = " << j1_refinement << ", j2_refinement = " << j2_refinement << std::endl;
    }
    j_refinement = j2_refinement;
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstructOnlyMultiScaling
                               (const CoefficientsByLevel<T> &u_scaling, int j,
                                CoefficientsByLevel<T> &u_loc_single, int &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    long k_refinement_first = 0L;
    for (const_coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        refCoeffs = basis.mra.phi.getRefinement(j,(*it).first,j_refinement,k_refinement_first);
        for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
            u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * (*it).second;
        }
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstructBSpline(int j, long k, T coeff,
                                                  CoefficientsByLevel<T> &u_loc_single,
                                                  int &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    j_refinement = 0;
    long k_refinement_first = 0L;
    refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,j_refinement,k_refinement_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstructWavelet(int j, long k, T coeff,
                                                  CoefficientsByLevel<T> &u_loc_single,
                                                  int &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    j_refinement = 0;
    long k_refinement_first = 0L;
    refCoeffs = basis.psi.getRefinement(j,k,j_refinement,k_refinement_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}


template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_(const CoefficientsByLevel<T>  &u_loc_single,
                                          CoefficientsByLevel<T>  &u_bspline, int j_bspline,
                                          CoefficientsByLevel<T>  &u_wavelet, int j_wavelet) const
{
    if (u_loc_single.map.size()==0) return;
    for (coeffbylevel_it it=u_bspline.map.begin(); it!=u_bspline.map.end(); ++it) {
        T coeff = this->decompose_BSpline(u_loc_single, j_bspline, (*it).first);
        (*it).second += coeff;
    }
    for (coeffbylevel_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        T coeff = this->decompose_Wavelet(u_loc_single, j_wavelet, (*it).first);
        (*it).second += coeff;
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_OnlyMultiScaling(const CoefficientsByLevel<T>  &u_loc_single,
                                                          CoefficientsByLevel<T>  &u_scaling, int j_scaling)
                                                          const
{
    if (u_loc_single.map.size()==0) return;
    for (coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        T coeff = this->decompose_Scaling(u_loc_single, j_scaling, (*it).first);
        (*it).second += coeff;
    }
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_Scaling(const CoefficientsByLevel<T> &u_loc_single,
                                                 int j, long k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    int refinement_j = 0;
    long refinement_k_first = 0L;
    T val = 0.;
    refCoeffs = basis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }
    return val;
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_BSpline(const CoefficientsByLevel<T> &u_loc_single,
                                                 int j, long k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    int refinement_j = 0;
    long refinement_k_first = 0L;
    T val = 0.;
    refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }
    return val;
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_Wavelet(const CoefficientsByLevel<T> &u_loc_single,
                                                 int j, long k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    int refinement_j = 0;
    long refinement_k_first = 0L;
    T val = 0.;
    refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }
    return val;
}

}   // namespace lawa
