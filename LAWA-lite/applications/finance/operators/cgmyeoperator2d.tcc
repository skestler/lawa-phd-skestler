namespace lawa {

template <typename Basis2D>
FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>::FinanceOperator2D
(const Basis2D& _basis, const ProcessParameters2D<T,CGMYeUnivariateJump2D> &_processparameters,
 T _R1_1, T _R2_1, T _R1_2, T _R2_2, int order,
 int _internal_compression_level1, int _internal_compression_level2)
: basis(_basis), processparameters(_processparameters),
  R1_1(_R1_1), R2_1(_R2_1), R1_2(_R1_2), R2_2(_R2_2),
  internal_compression_level1(_internal_compression_level1),
  internal_compression_level2(_internal_compression_level2),
  cgmyeop1d_1(basis.first, processparameters.proc_param1,R1_1,R2_1,order,internal_compression_level1),
  cgmyeop1d_2(basis.second,processparameters.proc_param2,R1_2,R2_2,order,internal_compression_level2),
  integral1(basis.first,basis.first), integral2(basis.second,basis.second)
{

}

template <typename Basis2D>
void
FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>::setCompressionLevel(int _internal_compression_level1,
                                                                       int _internal_compression_level2)
{
    internal_compression_level1 = _internal_compression_level1;
    cgmyeop1d_1.setCompressionLevel(internal_compression_level1);
    internal_compression_level2 = _internal_compression_level2;
    cgmyeop1d_2.setCompressionLevel(internal_compression_level2);
}

template <typename Basis2D>
typename Basis2D::T
FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>::operator()(const Index2D &row,
                                                              const Index2D &col)
{
    T ret = 0.;
    // evaluate x1-direction for multiwavelet discretization -> if x2 indices are not identical, the
    // corresponding entry vanishes!!
    if (row.index2.j==col.index2.j && row.index2.xtype==col.index2.xtype && row.index2.k==col.index2.k) {
        ret += cgmyeop1d_1(row.index1,col.index1);
    }
    // evaluate x2-direction for multiwavelet discretization -> if x1 indices are not identical, the
    // corresponding entry vanishes!!
    if (row.index1.j==col.index1.j && row.index1.xtype==col.index1.xtype && row.index1.k==col.index1.k) {
        ret += cgmyeop1d_2(row.index2,col.index2);
    }
    return ret;
}

template <typename Basis2D>
void
FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>::eval
(Coefficients<Lexicographical,T,Index2D> &v, Coefficients<Lexicographical,T,Index2D> &Av,
 const char* evalType)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator            const_coeff2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::iterator                  coeff2d_it;
    typedef typename std::map<Index2D,int,lt<Lexicographical,Index2D> >::const_iterator const_indicesToInt_it;

    bool current_v_Is_v = true;
    if (strcmp(evalType,"galerkin")==0) {
        //std::cerr << "   FinanceOperator2D: Treating case Galerkin index sets." << std::endl;
        if (current_v.size()==v.size()) {
            for (const_coeff2d_it it=current_v.begin(); it!=current_v.end(); ++it) {
                if (v.find((*it).first)==v.end()) current_v_Is_v = false;
            }
        }
        else {
            current_v_Is_v = false;
        }

        if (current_v_Is_v) {
            //std::cerr << "   FinanceOperator2D: Treating case where sparse matrix is pre-computed." << std::endl;
            this->performMV(v, Av);
        }
        else {
            //std::cerr << "   FinanceOperator2D: Treating case where no sparse matrix is pre-computed." << std::endl;
            int N = v.size();
            stiffnessMatrix.resize(N,N);
            indicesToInt.clear();
            current_v.clear();
            int row_count = 1;
            for (const_coeff2d_it it_row=v.begin(); it_row!=v.end(); ++it_row) {
                indicesToInt[(*it_row).first] = row_count;
                ++row_count;
                current_v[(*it_row).first] = 0.;
            }
            row_count = 1;
            for (const_coeff2d_it it_row=v.begin(); it_row!=v.end(); ++it_row) {
                int col_count = 1;
                for (const_coeff2d_it it_col=v.begin(); it_col!=v.end(); ++it_col) {
                    T val = this->operator()((*it_row).first,(*it_col).first);
                    if (fabs(val)>0) stiffnessMatrix(row_count,col_count) = val;
                    ++col_count;
                }
                ++row_count;
            }
            stiffnessMatrix.finalize();
            //std::stringstream matrixfilename;
            //matrixfilename << "_A_" << v.size() << "_" << processparameters;
            //spy(stiffnessMatrix,matrixfilename.str().c_str(),true,(T)0.);
            //DenseMatrixT DenseA(N,N);
            //densify(cxxblas::NoTrans,stiffnessMatrix,DenseA);
            //std::cerr << "  A = " << DenseA << std::endl;
            this->performMV(v, Av);
        }
    }
    else {
        for (coeff2d_it it_row=Av.begin(); it_row!=Av.end(); ++it_row) {
            T row_val = 0.;
            for (const_coeff2d_it it_col=v.begin(); it_col!=v.end(); ++it_col) {
                T matrixentry = this->operator()((*it_row).first,(*it_col).first);
                row_val += matrixentry * (*it_col).second;
            }
            (*it_row).second += row_val;
        }
    }
}

template <typename Basis2D>
void
FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>::performMV
(Coefficients<Lexicographical,T,Index2D> &v, Coefficients<Lexicographical,T,Index2D> &Av)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator            const_coeff2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::iterator                  coeff2d_it;
    typedef typename std::map<Index2D,int,lt<Lexicographical,Index2D> >::const_iterator const_indicesToInt_it;


    if (v.size()!=Av.size()) {
        std::cerr << "FinanceOperator2D: Dimensions do not match." << std::endl; exit(1);
    }
    int N = v.size();
    DenseVectorT x(N), Ax(N);
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        const_indicesToInt_it pos_it = indicesToInt.find((*it).first);
        int pos = (*pos_it).second;
        x(pos) = (*it).second;
    }
    Ax = stiffnessMatrix * x;
    for (coeff2d_it it=Av.begin(); it!=Av.end(); ++it) {
        const_indicesToInt_it pos_it = indicesToInt.find((*it).first);
        int pos = (*pos_it).second;
        (*it).second += Ax(pos);
    }

}

}   // namespace lawa
