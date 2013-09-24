namespace lawa {

template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
SparseGrid2D<T, IndexOneD, Basis2D, BilinearForm, RHS>::SparseGrid2D(const Basis2D &_basis,
                                                        BilinearForm &_a, RHS &_rhs, T _J, T _gamma)
: basis(_basis), a(_a), rhs(_rhs),
  J(_J), gamma(_gamma), Lambda(basis.first.d,basis.first.d_)
{
    setupIndexSet();
}

template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
IndexSet<Index2D>
SparseGrid2D<T, IndexOneD, Basis2D, BilinearForm, RHS>::getIndexSet()
{
    return Lambda;
}

template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
void
SparseGrid2D<T, IndexOneD, Basis2D, BilinearForm, RHS>::solve_cg(T &energynorm)
{
    int N = Lambda.size();
    DenseVector<Array<T> > b(N), Au(N);
    u.engine().resize(N);
    int row_count=1;
    std::cerr << "Right-hand side assembling started." << std::endl;
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
        b(row_count) = rhs(*row);
    }
    std::cerr << "Right-hand side assembling finished." << std::endl;
    flens::SparseGeMatrix<CRS<T,CRS_General> > A(N,N);
    lawa::toFlensSparseMatrix(a, Lambda, Lambda, A);
    std::cerr << "Sparse matrix assembled: (" << A.numNonZeros()
              << " / " << Lambda.size()*Lambda.size() << ")"<< std::endl;
    int NumOfIterations = lawa::cg(A,u,b,1e-14,1000);
    Au = A*u;
    energynorm = std::sqrt(u*Au);
    std::cerr << "SparseGrid2D::solve_cg required " << NumOfIterations << std::endl;
}


template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
T
SparseGrid2D<T, IndexOneD, Basis2D, BilinearForm, RHS>::evaluate(T x, T y)
{
    assert(int(u.length())==int(Lambda.size()));
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;

    T tmp=0.;
    int count=1;
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it,++count) {
        int j_x = (*it).index1.j,  k_x = (*it).index1.k;
        int j_y = (*it).index2.j,  k_y = (*it).index2.k;
        if (((*it).index1.xtype==XBSpline) && ((*it).index2.xtype==XBSpline)) {
            tmp+=a.prec(*it)*u(count)*basis.first.mra.phi(x,j_x,k_x,0)*basis.second.mra.phi(y,j_y,k_y,0);
        }
        else if (((*it).index1.xtype==XBSpline) && ((*it).index2.xtype==XWavelet)) {
            tmp+=a.prec(*it)*u(count)*basis.first.mra.phi(x,j_x,k_x,0)*basis.second.psi(y,j_y,k_y,0);
        }
        else if (((*it).index1.xtype==XWavelet) && ((*it).index2.xtype==XBSpline)) {
            tmp+=a.prec(*it)*u(count)*basis.first.psi(x,j_x,k_x,0)*basis.second.mra.phi(y,j_y,k_y,0);
        }
        else {
            tmp+=a.prec(*it)*u(count)*basis.first.psi(x,j_x,k_x,0)*basis.second.psi(y,j_y,k_y,0);
        }
    }
    return tmp;
}


template <typename T, typename IndexOneD, typename Basis2D, typename BilinearForm, typename RHS>
void
SparseGrid2D<T, IndexOneD, Basis2D, BilinearForm, RHS>::setupIndexSet()
{
    IndexSet<Index2D> ret(basis.first.d,basis.first.d_);
    int j0_x = basis.first.mra.j0;
    int j0_y = basis.second.mra.j0;

    for (int k_x=basis.first.mra.rangeI(j0_x).firstIndex();
             k_x<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k_x) {
        for (int k_y=basis.second.mra.rangeI(j0_y).firstIndex();
                 k_y<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k_y) {
            ret.insert(Index2D(IndexOneD(j0_x,k_x,XBSpline),IndexOneD(j0_y,k_y,XBSpline)));
        }
    }

    for (int k_x=basis.first.mra.rangeI(j0_x).firstIndex();
             k_x<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k_x) {
        for (int j_y=j0_y; j_y<=J+j0_y; ++j_y) {
            for (int k_y=basis.second.rangeJ(j_y).firstIndex();
                     k_y<=basis.second.rangeJ(j_y).lastIndex(); ++k_y) {
                ret.insert(Index2D(IndexOneD(j0_x,k_x,XBSpline),IndexOneD(j_y,k_y,XWavelet)));
            }
        }
    }

    for (int j_x=j0_x; j_x<=J+j0_x; ++j_x) {
        for (int k_x=basis.first.rangeJ(j_x).firstIndex();
                 k_x<=basis.first.rangeJ(j_x).lastIndex(); ++k_x) {
            for (int k_y=basis.second.mra.rangeI(j0_y).firstIndex();
                     k_y<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k_y) {
                ret.insert(Index2D(IndexOneD(j_x,k_x,XWavelet),IndexOneD(j0_y,k_y,XBSpline)));
            }
        }
    }

    for (int j_x=j0_x; j_x<=J+j0_x; ++j_x) {
        for (int j_y=j0_y; j_y<=J+j0_y; ++j_y) {
            if ((j_x-j0_x)+(j_y-j0_y) > (1-gamma)*J+gamma*std::max(j_x-j0_x,j_y-j0_y)) continue;

            for (int k_x=basis.first.rangeJ(j_x).firstIndex();
                     k_x<=basis.first.rangeJ(j_x).lastIndex(); ++k_x) {
                for (int k_y=basis.second.rangeJ(j_y).firstIndex();
                         k_y<=basis.second.rangeJ(j_y).lastIndex(); ++k_y) {
                    ret.insert(Index2D(IndexOneD(j_x,k_x,XWavelet),IndexOneD(j_y,k_y,XWavelet)));
                }
            }
        }
    }
    Lambda=ret;
}


}   //namespace lawa
