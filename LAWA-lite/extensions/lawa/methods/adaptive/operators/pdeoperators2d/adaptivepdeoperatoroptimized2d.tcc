namespace lawa {

template <typename T, DomainType Domain1, DomainType Domain2>
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::AdaptivePDEOperatorOptimized2D(const Basis2D &_basis, T _reaction, T _convection_x,
                                 T _convection_y, T _diffusion_y)
: basis(_basis), reaction(_reaction), convection_x(_convection_x), convection_y(_convection_y),
  diffusion_y(_diffusion_y),
  cA(0.), CA(0.), kappa(0.),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  laplace_data1d_x(basis.first), convection_data1d_x(basis.first), identity_data1d_x(basis.first),
  laplace_data1d_y(basis.second), convection_data1d_y(basis.second), identity_data1d_y(basis.second),
  P_data()
{
    T cA_x=0., CA_x=0., cA_y=0., CA_y = 0.;
    std::cout << "AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,"
              << "Primal,Domain2,SparseMulti>: j0_x=" << basis.first.j0 << ", j0_y="
              << basis.second.j0 << ", "
              << "reaction=" << reaction << ", convection_x=" << convection_x << ", "
              << "convection_y=" << convection_y << std::endl;

    if (basis.first.d==4 && basis.second.d==4) {

        if      (basis.first.j0==0)  {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0== 1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-2) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-3) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-4) {    cA_x = 0.19;  CA_x = 2.9;    }
        else assert(0);

        if      (basis.second.j0==0)  {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0== 1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-2) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-3) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-4) {    cA_y = 0.19;  CA_y = 2.9;    }
        else assert(0);

    }
    cA = std::min(cA_x*cA_x,cA_y*cA_y);
    CA = std::max(CA_x*CA_x,CA_y*CA_y);
    kappa = CA/cA;
}

template <typename T, DomainType Domain1, DomainType Domain2>
T
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::operator()(const Index2D &row_index, const Index2D &col_index)
{
    T id_x = identity_data1d_x(row_index.index1,col_index.index1);
    T id_y = identity_data1d_y(row_index.index2,col_index.index2);

    T dd_y=0., d_y=0.;
    if (fabs(id_x)>0) {
        dd_y = laplace_data1d_y(row_index.index2,col_index.index2);
        d_y  = convection_data1d_y(row_index.index2,col_index.index2);
    }
    T dd_x=0., d_x=0.;
    if (fabs(id_y)>0) {
        dd_x = laplace_data1d_x(row_index.index1,col_index.index1);
        d_x  = convection_data1d_x(row_index.index1,col_index.index1);
    }

    T val =   (            dd_x + convection_x*d_x + 0.5*reaction*id_x)*id_y
            + (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y)*id_x;

    val *= this->prec(row_index) * this->prec(col_index);
    return val;
}

template <typename T, DomainType Domain1, DomainType Domain2>
T
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::prec(const Index2D &index)
{
    T prec = 1.;
/*
    const_coeff2d_it it_P_end   = P_data.end();
    const_coeff2d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
*/
        T prec_id_x = identity_data1d_x(index.index1,index.index1);
        T prec_id_y = identity_data1d_y(index.index2,index.index2);
        T prec_dd_x = laplace_data1d_x(index.index1,index.index1);
        T prec_dd_y = laplace_data1d_y(index.index2,index.index2);
        T tmp = 1./std::sqrt(fabs(   (            prec_dd_x + 0.5*reaction*prec_id_x)*prec_id_y
                                   + (diffusion_y*prec_dd_y + 0.5*reaction*prec_id_y)*prec_id_x ) );
//        P_data[index] = tmp;
        prec *= tmp;
//    }
    return prec;
}

template <typename T, DomainType Domain1, DomainType Domain2>
Coefficients<Lexicographical,T,Index2D>
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &v)
{
    Coefficients<Lexicographical,T,Index2D> ret;

    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index2D col_index = (*col).first;
        T prec_val_col_index = this->prec(col_index) * (*col).second;

        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        if (sparsitypatterns_x.count(col_index.index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern(col_index.index1, LambdaRow_x);
            sparsitypatterns_x[col_index.index1] = LambdaRowSparse_x;
        }
        else {
            LambdaRowSparse_x = sparsitypatterns_x[col_index.index1];
        }

        if (sparsitypatterns_y.count(col_index.index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern(col_index.index2, LambdaRow_y);
            sparsitypatterns_y[col_index.index2] = LambdaRowSparse_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[col_index.index2];
        }

        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            T id_x = identity_data1d_x(*row_x,col_index.index1);
            T d_x  = convection_data1d_x(*row_x,col_index.index1);
            T dd_x = laplace_data1d_x(*row_x,col_index.index1);
            if (fabs(id_x)<1e-13 && fabs(d_x)<1e-12 && fabs(dd_x)<1e-13) continue;
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            Index2D row_index(*row_x,*row_y);
                if (LambdaRow.count(row_index)>0) {
                    T id_y = identity_data1d_y(row_index.index2,col_index.index2);
                    T d_y  = convection_data1d_y(row_index.index2,col_index.index2);
                    T dd_y = laplace_data1d_y(row_index.index2,col_index.index2);


                    T val =   (            dd_x + convection_x*d_x + 0.5*reaction*id_x)*id_y
                            + (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y)*id_x;

                    if (fabs(val)>0.) {
                        ret[row_index] += val * prec_val_col_index;
                    }
                }
            }
        }

    }
    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *=  this->prec((*it).first);
    }
    return ret;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, int J)
{
    std::cerr << "  -> toFlensSparseMatrix called with J= " << J << std::endl;

    int maxSizeSparsityPattern=0;

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;
    std::map<Index1D,Coefficients<Lexicographical,T,Index1D>,lt<Lexicographical,Index1D> >
        sparsitypatterns_identity_x, sparsitypatterns_convection_x, sparsitypatterns_laplace_x,
        sparsitypatterns_identity_y, sparsitypatterns_convection_y, sparsitypatterns_laplace_y;

    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);


    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;

    int row_count = 1;
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    std::map<Index2D,int,lt<Lexicographical,Index2D> >::const_iterator row_indices_end;
    row_indices_end = row_indices.end();
    const_set2d_it LambdaRow_end=LambdaRow.end();

    int col_count = 1;
    for (const_set2d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        Index2D col_index = *col;
        T prec_col_index = this->prec(col_index);

        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        Coefficients<Lexicographical,T,Index1D>
            LambdaRowSparseIdentity_x, LambdaRowSparseConvection_x, LambdaRowSparseLaplace_x,
            LambdaRowSparseIdentity_y, LambdaRowSparseConvection_y, LambdaRowSparseLaplace_y;
        if (sparsitypatterns_x.count((*col).index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern((*col).index1, LambdaRow_x, 1);
            sparsitypatterns_x[(*col).index1] = LambdaRowSparse_x;
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                LambdaRowSparseIdentity_x[*row_x]   = identity_data1d_x(*row_x,col_index.index1);
                LambdaRowSparseConvection_x[*row_x] = convection_data1d_x(*row_x,col_index.index1);
                LambdaRowSparseLaplace_x[*row_x]    = laplace_data1d_x(*row_x,col_index.index1);
            }
            sparsitypatterns_identity_x[(*col).index1]   = LambdaRowSparseIdentity_x;
            sparsitypatterns_convection_x[(*col).index1] = LambdaRowSparseConvection_x;
            sparsitypatterns_laplace_x[(*col).index1]    = LambdaRowSparseLaplace_x;
        }
        else {
            LambdaRowSparse_x           = sparsitypatterns_x[(*col).index1];
            LambdaRowSparseIdentity_x   = sparsitypatterns_identity_x[(*col).index1];
            LambdaRowSparseConvection_x = sparsitypatterns_convection_x[(*col).index1];
            LambdaRowSparseLaplace_x    = sparsitypatterns_laplace_x[(*col).index1];
        }

        if (sparsitypatterns_y.count((*col).index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern((*col).index2, LambdaRow_y, 1);
            sparsitypatterns_y[(*col).index2] = LambdaRowSparse_y;
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                LambdaRowSparseIdentity_y[*row_y]   = identity_data1d_y(*row_y,col_index.index2);
                LambdaRowSparseConvection_y[*row_y] = convection_data1d_y(*row_y,col_index.index2);
                LambdaRowSparseLaplace_y[*row_y]    = laplace_data1d_y(*row_y,col_index.index2);
            }
            sparsitypatterns_identity_y[(*col).index2] = LambdaRowSparseIdentity_y;
            sparsitypatterns_convection_y[(*col).index2] = LambdaRowSparseConvection_y;
            sparsitypatterns_laplace_y[(*col).index2]  = LambdaRowSparseLaplace_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[(*col).index2];
            LambdaRowSparseIdentity_y   = sparsitypatterns_identity_y[(*col).index2];
            LambdaRowSparseConvection_y = sparsitypatterns_convection_y[(*col).index2];
            LambdaRowSparseLaplace_y    = sparsitypatterns_laplace_y[(*col).index2];
        }

        maxSizeSparsityPattern=std::max(maxSizeSparsityPattern,(int)(LambdaRowSparse_x.size()*LambdaRowSparse_y.size()));


        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            T id_x = LambdaRowSparseIdentity_x[*row_x];
            T d_x  = LambdaRowSparseConvection_x[*row_x];
            T dd_x = LambdaRowSparseLaplace_x[*row_x];
            if (fabs(id_x)<1e-13 && fabs(d_x)<1e-13 && fabs(dd_x)<1e-13) continue;
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                Index2D row_index(*row_x,*row_y);
                std::map<Index2D,int,lt<Lexicographical,Index2D> >::const_iterator row_count_ptr;
                row_count_ptr = row_indices.find(row_index);
                if (row_count_ptr!=row_indices_end) {
                    T id_y = LambdaRowSparseIdentity_y[*row_y];
                    T d_y  = LambdaRowSparseConvection_y[*row_y];
                    T dd_y = LambdaRowSparseLaplace_y[*row_y];

                    T val =   (            dd_x + convection_x*d_x + 0.5*reaction*id_x)*id_y
                            + (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y)*id_x;

                    T prec_row_index = this->prec(row_index);
                    T prec_val = prec_row_index* val * prec_col_index;
                    if (fabs(prec_val)>0.) {
                        A_flens((*row_count_ptr).second,col_count) = prec_val;
                    }
                }
            }
        }

    }
    std::cout << "Size of LambdaRow: " << LambdaRow.size()
              << ", size of LambdaRowSparse: " << maxSizeSparsityPattern << std::endl;

    A_flens.finalize();
    std::cerr << "   Number of non-zeros: " << A_flens.numNonZeros() << std::endl;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, T /*eps*/)
{
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,1);
}

template <typename T, DomainType Domain1, DomainType Domain2>
Coefficients<Lexicographical,T,Index2D>
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, int /*k*/, int /*J*/,
        cxxblas::Transpose trans)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    this->apply(v,0.,ret,trans);
    return ret;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
        Coefficients<Lexicographical,T,Index2D> &ret, cxxblas::Transpose trans)
{
    if (v.size()==0) return;

    Index1D_Coefficients1D_Hash y_v;
    Coefficients<Lexicographical,T,Index2D> I_S_v(SIZEHASHINDEX2D);
    Index1D_Coefficients1D_Hash x_I_S_v;

    Timer time;
    time.start();

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        if (y_v.count((*col).first.index2)>0) {
            y_v[(*col).first.index2].operator[]((*col).first.index1) = this->prec((*col).first) * (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_x;
            coeff_x[(*col).first.index1] = this->prec((*col).first) * (*col).second;
            y_v[(*col).first.index2] = coeff_x;
        }
    }

    if (trans==cxxblas::NoTrans) {

        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D col_index_y = (*it).first;

            IndexSet<Index1D> LambdaRowSparse_y;
            LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                              col_index_y.j+1, false);
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                T dd_y = laplace_data1d_y(*row_y,col_index_y);
                T id_y = identity_data1d_y(*row_y,col_index_y);
                T d_y  = convection_data1d_y(*row_y,col_index_y);
                T val2_y = (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y);

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index2D row_index((*coeff_x_it).first,*row_y);
                    T val2 = val2_y * (*coeff_x_it).second;
                    if (fabs(val2)>0) I_S_v[row_index] += val2;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }
        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            Index1D col_index_x = (*it).first;
            IndexSet<Index1D> LambdaRowSparse_x;
            LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                                  col_index_x.j+1, false);
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                T id_x = identity_data1d_x(*row_x,col_index_x);
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index2D row_index(*row_x,(*coeff_y_it).first);
                    T val = id_x * (*coeff_y_it).second;
                    if (fabs(val)>0.) ret[row_index] += val;
                }
            }
        }

        for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
                    (*it).second = 0.;
        }
        for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
                (*coeff_it).second = 0.;
            }
        }


        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D col_index_y = (*it).first;

            IndexSet<Index1D> LambdaRowSparse_y;
            LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                              col_index_y.j+1, false);
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                T id_y = identity_data1d_y(*row_y,col_index_y);
                T val1_y =  id_y;

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index2D row_index((*coeff_x_it).first,*row_y);
                    T val1 = val1_y * (*coeff_x_it).second;
                    if (fabs(val1)>0) I_S_v[row_index] += val1;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }

        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

            Index1D col_index_x = (*it).first;
            IndexSet<Index1D> LambdaRowSparse_x;
            LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                                  col_index_x.j+1, false);
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                T id_x = identity_data1d_x(*row_x,col_index_x);
                T dd_x = laplace_data1d_x(*row_x,col_index_x);
                T d_x  = convection_data1d_x(*row_x,col_index_x);
                T val_x = dd_x + convection_x*d_x + 0.5*reaction*id_x;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index2D row_index(*row_x,(*coeff_y_it).first);
                    T val = val_x * (*coeff_y_it).second;
                    if (fabs(val)>0.) ret[row_index] += val;
                }
            }
        }


        for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {

            (*it).second *=  this->prec((*it).first);
        }

        time.stop();
        //std::cerr << "      New structure required: " << time.elapsed() << std::endl;

        //std::cerr << "      Input length = " << v.size()  << ", output length = " << ret.size()
        //          << ", quotient(output vs. input) = " << T(ret.size())/T(v.size()) << std::endl;
    }
    else if (trans==cxxblas::Trans) {

        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D row_index_y = (*it).first;

            IndexSet<Index1D> LambdaColSparse_y;

            LambdaColSparse_y = lambdaTilde1d_PDE(row_index_y, basis.second, 1, std::max(row_index_y.j-1,basis.second.j0),
                                                              row_index_y.j+1, false);
            for (const_set1d_it col_y=LambdaColSparse_y.begin(); col_y!=LambdaColSparse_y.end(); ++col_y) {
                T dd_y = laplace_data1d_y(row_index_y,*col_y);
                T id_y = identity_data1d_y(row_index_y,*col_y);
                T d_y  = convection_data1d_y(row_index_y,*col_y);
                T val2_y = (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y);

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index1D row_index_x = (*coeff_x_it).first;
                    Index2D col_index(row_index_x,*col_y);
                    T val2 = val2_y * (*coeff_x_it).second;
                    if (fabs(val2)>0) I_S_v[col_index] += val2;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }
        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            Index1D row_index_x = (*it).first;
            IndexSet<Index1D> LambdaColSparse_x;
            LambdaColSparse_x = lambdaTilde1d_PDE(row_index_x, basis.first, 1, std::max(row_index_x.j-1,basis.first.j0),
                                                  row_index_x.j+1, false);
            for (const_set1d_it col_x=LambdaColSparse_x.begin(); col_x!=LambdaColSparse_x.end(); ++col_x) {
                T id_x = identity_data1d_x(row_index_x,*col_x);
                if (!(fabs(id_x)>0)) continue;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index1D row_index_y = (*coeff_y_it).first;
                    Index2D col_index(*col_x,row_index_y);
                    T val = id_x * (*coeff_y_it).second;
                    ret[col_index] += val;
                }
            }
        }

        for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
            (*it).second = 0.;
        }
        for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
                (*coeff_it).second = 0.;
            }
        }


        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D row_index_y = (*it).first;

            IndexSet<Index1D> LambdaColSparse_y;

            LambdaColSparse_y = lambdaTilde1d_PDE(row_index_y, basis.second, 1, std::max(row_index_y.j-1,basis.second.j0),
                                                              row_index_y.j+1, false);
            for (const_set1d_it col_y=LambdaColSparse_y.begin(); col_y!=LambdaColSparse_y.end(); ++col_y) {
                T id_y = identity_data1d_y(row_index_y,*col_y);
                T val1_y =  id_y;

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index1D row_index_x = (*coeff_x_it).first;
                    Index2D col_index(row_index_x,*col_y);
                    T val1 = val1_y * (*coeff_x_it).second;
                    if (fabs(val1)>0) I_S_v[col_index] += val1;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }

        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

            Index1D row_index_x = (*it).first;
            IndexSet<Index1D> LambdaColSparse_x;
            Coefficients<Lexicographical,T,Index1D> Colvalues_A_x_index;
            LambdaColSparse_x = lambdaTilde1d_PDE(row_index_x, basis.first, 1, std::max(row_index_x.j-1,basis.first.j0),
                                                  row_index_x.j+1, false);
            for (const_set1d_it col_x=LambdaColSparse_x.begin(); col_x!=LambdaColSparse_x.end(); ++col_x) {
                T id_x = identity_data1d_x(row_index_x,*col_x);
                T dd_x = laplace_data1d_x(row_index_x,*col_x);
                T d_x  = convection_data1d_x(row_index_x,*col_x);
                T val_x = dd_x + convection_x*d_x + 0.5*reaction*id_x;
                if (!(fabs(val_x)>0)) continue;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index1D row_index_y = (*coeff_y_it).first;
                    Index2D col_index(*col_x,row_index_y);
                    T val = val_x * (*coeff_y_it).second;
                    ret[col_index] += val;
                }
            }
        }

        for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
            (*it).second *=  this->prec((*it).first);
        }

        time.stop();
        //std::cerr << "      A^T (A v) required: " << time.elapsed() << std::endl;
        //std::cerr << "      Size of P_data: " << P_data.size() << std::endl;
        //std::cerr << "      Input length = " << v.size()  << ", output length = " << ret.size()
        //          << ", quotient(output vs. input) = " << T(ret.size())/T(v.size()) << std::endl;

    }

    y_v.clear();
    I_S_v.clear();
    for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        (*it).second.clear();
    }
    x_I_S_v.clear();


}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
        const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
        cxxblas::Transpose trans)
{
    if (v.size()==0) return;

    IndexSet<Index1D> Lambda_x, Lambda_y;
    split(Lambda, Lambda_x, Lambda_y);

    Index1D_Coefficients1D_Hash y_v;
    Coefficients<Lexicographical,T,Index2D> I_S_v(SIZEHASHINDEX2D);
    Index1D_Coefficients1D_Hash x_I_S_v;

    Timer time;
    time.start();

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        if (y_v.count((*col).first.index2)>0) {
            y_v[(*col).first.index2].operator[]((*col).first.index1) = this->prec((*col).first) * (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_x;
            coeff_x[(*col).first.index1] = this->prec((*col).first) * (*col).second;
            y_v[(*col).first.index2] = coeff_x;
        }
    }

    if (trans==cxxblas::NoTrans) {

        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D col_index_y = (*it).first;

            IndexSet<Index1D> LambdaRowSparse_y;
            LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                              col_index_y.j+1, false);
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                if (Lambda_y.count(*row_y)==0) continue;
                T dd_y = laplace_data1d_y(*row_y,col_index_y);
                T id_y = identity_data1d_y(*row_y,col_index_y);
                T d_y  = convection_data1d_y(*row_y,col_index_y);
                T val2_y = (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y);

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index2D row_index((*coeff_x_it).first,*row_y);
                    T val2 = val2_y * (*coeff_x_it).second;
                    if (fabs(val2)>0) I_S_v[row_index] += val2;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }
        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            Index1D col_index_x = (*it).first;
            IndexSet<Index1D> LambdaRowSparse_x;
            LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                                  col_index_x.j+1, false);
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                if (Lambda_x.count(*row_x)==0) continue;
                T id_x = identity_data1d_x(*row_x,col_index_x);
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index2D row_index(*row_x,(*coeff_y_it).first);
                    if (Lambda.count(row_index)==0) continue;
                    T val = id_x * (*coeff_y_it).second;
                    if (fabs(val)>0.) ret[row_index] += val;
                }
            }
        }

        for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
                    (*it).second = 0.;
        }
        for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
                (*coeff_it).second = 0.;
            }
        }


        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D col_index_y = (*it).first;

            IndexSet<Index1D> LambdaRowSparse_y;
            LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                              col_index_y.j+1, false);
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                if (Lambda_y.count(*row_y)==0) continue;
                T id_y = identity_data1d_y(*row_y,col_index_y);
                T val1_y =  id_y;

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index2D row_index((*coeff_x_it).first,*row_y);
                    T val1 = val1_y * (*coeff_x_it).second;
                    if (fabs(val1)>0) I_S_v[row_index] += val1;
                }
            }
        }

        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }

        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

            Index1D col_index_x = (*it).first;
            IndexSet<Index1D> LambdaRowSparse_x;
            LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                                  col_index_x.j+1, false);
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                if (Lambda_x.count(*row_x)==0) continue;
                T id_x = identity_data1d_x(*row_x,col_index_x);
                T dd_x = laplace_data1d_x(*row_x,col_index_x);
                T d_x  = convection_data1d_x(*row_x,col_index_x);
                T val_x = dd_x + convection_x*d_x + 0.5*reaction*id_x;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index2D row_index(*row_x,(*coeff_y_it).first);
                    if (Lambda.count(row_index)==0) continue;
                    T val = val_x * (*coeff_y_it).second;
                    if (fabs(val)>0.) ret[row_index] += val;
                }
            }
        }
        for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {

            (*it).second *=  this->prec((*it).first);
        }
        time.stop();
        //std::cerr << "      New structure required: " << time.elapsed() << std::endl;

        //std::cerr << "      Input length = " << v.size()  << ", output length = " << ret.size()
        //          << ", quotient(output vs. input) = " << T(ret.size())/T(v.size()) << std::endl;
    }
    else if (trans==cxxblas::Trans) {

        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D row_index_y = (*it).first;

            IndexSet<Index1D> LambdaColSparse_y;

            LambdaColSparse_y = lambdaTilde1d_PDE(row_index_y, basis.second, 1, std::max(row_index_y.j-1,basis.second.j0),
                                                              row_index_y.j+1, false);
            for (const_set1d_it col_y=LambdaColSparse_y.begin(); col_y!=LambdaColSparse_y.end(); ++col_y) {
                if (Lambda_y.count(*col_y)==0) continue;
                T dd_y = laplace_data1d_y(row_index_y,*col_y);
                T id_y = identity_data1d_y(row_index_y,*col_y);
                T d_y  = convection_data1d_y(row_index_y,*col_y);
                T val2_y = (diffusion_y*dd_y + convection_y*d_y + 0.5*reaction*id_y);

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index1D row_index_x = (*coeff_x_it).first;
                    Index2D col_index(row_index_x,*col_y);
                    T val2 = val2_y * (*coeff_x_it).second;
                    if (fabs(val2)>0) I_S_v[col_index] += val2;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }
        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            Index1D row_index_x = (*it).first;
            IndexSet<Index1D> LambdaColSparse_x;
            LambdaColSparse_x = lambdaTilde1d_PDE(row_index_x, basis.first, 1, std::max(row_index_x.j-1,basis.first.j0),
                                                  row_index_x.j+1, false);
            for (const_set1d_it col_x=LambdaColSparse_x.begin(); col_x!=LambdaColSparse_x.end(); ++col_x) {
                if (Lambda_x.count(*col_x)==0) continue;
                T id_x = identity_data1d_x(row_index_x,*col_x);
                if (!(fabs(id_x)>0)) continue;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index1D row_index_y = (*coeff_y_it).first;
                    Index2D col_index(*col_x,row_index_y);
                    if (Lambda.count(col_index)==0) continue;
                    T val = id_x * (*coeff_y_it).second;
                    ret[col_index] += val;
                }
            }
        }

        for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
            (*it).second = 0.;
        }
        for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
            for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
                (*coeff_it).second = 0.;
            }
        }


        for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
            Index1D row_index_y = (*it).first;

            IndexSet<Index1D> LambdaColSparse_y;

            LambdaColSparse_y = lambdaTilde1d_PDE(row_index_y, basis.second, 1, std::max(row_index_y.j-1,basis.second.j0),
                                                              row_index_y.j+1, false);
            for (const_set1d_it col_y=LambdaColSparse_y.begin(); col_y!=LambdaColSparse_y.end(); ++col_y) {
                if (Lambda_y.count(*col_y)==0) continue;
                T id_y = identity_data1d_y(row_index_y,*col_y);
                T val1_y =  id_y;

                for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                    Index1D row_index_x = (*coeff_x_it).first;
                    Index2D col_index(row_index_x,*col_y);
                    T val1 = val1_y * (*coeff_x_it).second;
                    if (fabs(val1)>0) I_S_v[col_index] += val1;
                }
            }
        }
        for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
            if (x_I_S_v.count((*col).first.index1)>0) {
                x_I_S_v[(*col).first.index1].operator[]((*col).first.index2) = (*col).second;
            }
            else {
                Coefficients<Lexicographical,T,Index1D> coeff_y;
                coeff_y[(*col).first.index2] = (*col).second;
                x_I_S_v[(*col).first.index1] = coeff_y;
            }
        }
        for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

            Index1D row_index_x = (*it).first;
            IndexSet<Index1D> LambdaColSparse_x;
            Coefficients<Lexicographical,T,Index1D> Colvalues_A_x_index;
            LambdaColSparse_x = lambdaTilde1d_PDE(row_index_x, basis.first, 1, std::max(row_index_x.j-1,basis.first.j0),
                                                  row_index_x.j+1, false);
            for (const_set1d_it col_x=LambdaColSparse_x.begin(); col_x!=LambdaColSparse_x.end(); ++col_x) {
                if (Lambda_x.count(*col_x)==0) continue;
                T id_x = identity_data1d_x(row_index_x,*col_x);
                T dd_x = laplace_data1d_x(row_index_x,*col_x);
                T d_x  = convection_data1d_x(row_index_x,*col_x);
                T val_x = dd_x + convection_x*d_x + 0.5*reaction*id_x;
                if (!(fabs(val_x)>0)) continue;
                for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                    Index1D row_index_y = (*coeff_y_it).first;
                    Index2D col_index(*col_x,row_index_y);
                    if (Lambda.count(col_index)==0) continue;
                    T val = val_x * (*coeff_y_it).second;
                    ret[col_index] += val;
                }
            }
        }

        for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {

            (*it).second *=  this->prec((*it).first);
        }

        time.stop();
        //std::cerr << "      A^T (A v) required: " << time.elapsed() << std::endl;
        //std::cerr << "      Size of P_data: " << P_data.size() << std::endl;
        //std::cerr << "      Input length = " << v.size()  << ", output length = " << ret.size()
        //          << ", quotient(output vs. input) = " << T(ret.size())/T(v.size()) << std::endl;

    }

    y_v.clear();
    I_S_v.clear();
    for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        (*it).second.clear();
    }
    x_I_S_v.clear();

}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptivePDEOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::clear()
{
    Coefficients<Lexicographical,T,Index2D> tmp;
    P_data = tmp;
}

}   //lawa
