namespace lawa {

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index1D>  &v,
                Coefficients<Lexicographical,T,Index1D>  &C_v, const char* residualType,
                bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index;
        C_index = C((*it).first, (T)1., basis);
        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it newindex=C_index.begin(); newindex!=C_index.end(); ++newindex) {
                if (C_v.find(*newindex)==C_v.end()) {
                    completeMultiTree(basis,*newindex,C_v, sparsetree);
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }
    return;
}

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v, const char* residualType,
                bool IsMW, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index2D>::iterator                                      set2d_it;


#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif BOOST
    typedef boost::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer(v.size());

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2;

        //C_index1 = C((*it).first.index1, (T)1., basis.first);
        //C_index2 = C((*it).first.index2, (T)1., basis.second);

        IndexConeContainer::const_iterator indexConeContainer_ptr;
        indexConeContainer_ptr = indexConeContainer.find((*it).first.index1);
        if (indexConeContainer_ptr!=indexConeContainer.end()) {
            C_index1 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            indexConeContainer[(*it).first.index1] = C_index1;
        }
        indexConeContainer_ptr = indexConeContainer.find((*it).first.index2);
        if (indexConeContainer_ptr!=indexConeContainer.end()) {
            C_index2 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            indexConeContainer[(*it).first.index2] = C_index2;
        }

        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index2D newindex((*it_C_index1), (*it).first.index2);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,1,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D newindex((*it).first.index1, (*it_C_index2));
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,2,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
        else if (strcmp(residualType,"large1")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    Index2D newindex((*it_C_index1), (*it_C_index2));
                    if (C_v.find(newindex)==C_v.end()) {
                        completeMultiTree(basis,newindex,C_v,0,sparsetree);
                    }
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }


    return;
}

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                Coefficients<Lexicographical,T,Index3D>  &C_v, const char* residualType,
                bool IsMW, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index3D>::const_iterator                                set3d_it;

#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif BOOST
    typedef boost::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer(v.size());

    for (const_coeff3d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2, C_index3;
        IndexConeContainer::const_iterator indexConeContainer_ptr;
        indexConeContainer_ptr = indexConeContainer.find((*it).first.index1);
        if (indexConeContainer_ptr!=indexConeContainer.end()) {
            C_index1 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            indexConeContainer[(*it).first.index1] = C_index1;
        }
        indexConeContainer_ptr = indexConeContainer.find((*it).first.index2);
        if (indexConeContainer_ptr!=indexConeContainer.end()) {
            C_index2 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            indexConeContainer[(*it).first.index2] = C_index2;
        }
        indexConeContainer_ptr = indexConeContainer.find((*it).first.index3);
        if (indexConeContainer_ptr!=indexConeContainer.end()) {
            C_index3 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index3 = C((*it).first.index3, (T)1., basis.second);
            indexConeContainer[(*it).first.index3] = C_index3;
        }
        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index3D newindex((*it_C_index1), (*it).first.index2, (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,1,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index3D newindex((*it).first.index1, (*it_C_index2), (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,2,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
                Index3D newindex((*it).first.index1, (*it).first.index2, (*it_C_index3));
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,3,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
        else if (strcmp(residualType,"large1")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    for (const_set1d_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
                        Index3D newindex((*it_C_index1), (*it_C_index2), (*it_C_index3));
                        if (C_v.find(newindex)==C_v.end()) {
                            completeMultiTree(basis,newindex,C_v,0,sparsetree);
                        }
                    }
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }
    return;
}

/*
template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                Coefficients<Lexicographical,T,Index3D>  &C_v, int coordDirec)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index3D>::const_iterator                                set3d_it;

#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer(v.size());

    for (const_coeff3d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        if (coordDirec == 1) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index1);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index1, (T)1., basis.first);
                indexConeContainer[(*it).first.index1] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it_C_index), (*it).first.index2, (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else if (coordDirec == 2) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index2);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index2, (T)1., basis.second);
                indexConeContainer[(*it).first.index2] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it).first.index1, (*it_C_index), (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else if (coordDirec == 3) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index3);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index3, (T)1., basis.second);
                indexConeContainer[(*it).first.index3] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it).first.index1, (*it).first.index2, (*it_C_index));
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: non-admissible coordinate direction " << coordDirec
                      << std::endl;
            exit(1);
        }
    }
    return;
}
*/

template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                          Coefficients<Lexicographical,T,Index2D>  &C_v, int J, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    IndexSet<Index1D> LambdaB_x1;
    for (int k= basis.first.mra.rangeI(basis.first.j0).firstIndex();
             k<=basis.first.mra.rangeI(basis.first.j0).lastIndex(); ++k) {
        LambdaB_x1.insert(Index1D(basis.first.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.first.rangeJL(j).firstIndex(); k<=basis.first.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.first.rangeJR(j).firstIndex(); k<=basis.first.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x2;
    for (int k= basis.second.mra.rangeI(basis.second.j0).firstIndex();
             k<=basis.second.mra.rangeI(basis.second.j0).lastIndex(); ++k) {
        LambdaB_x2.insert(Index1D(basis.second.j0,k,XBSpline));
    }
    for (int j=basis.second.j0; j<=J; ++j) {
        for (int k=basis.second.rangeJL(j).firstIndex(); k<=basis.second.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.second.rangeJR(j).firstIndex(); k<=basis.second.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
    }

    for (const_set1d_it it_x1=LambdaB_x1.begin(); it_x1!=LambdaB_x1.end(); ++it_x1) {
        for (const_set1d_it it_x2=LambdaB_x2.begin(); it_x2!=LambdaB_x2.end(); ++it_x2) {
            Index2D newindex((*it_x1), (*it_x2));
            if (C_v.find(newindex)==C_v.end()) {
                completeMultiTree(basis,newindex,C_v,0,sparsetree);
            }
        }
    }
}

template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                          Coefficients<Lexicographical,T,Index3D>  &C_v, int J, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    IndexSet<Index1D> LambdaB_x1;
    for (int k= basis.first.mra.rangeI(basis.first.j0).firstIndex();
             k<=basis.first.mra.rangeI(basis.first.j0).lastIndex(); ++k) {
        LambdaB_x1.insert(Index1D(basis.first.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.first.rangeJL(j).firstIndex(); k<=basis.first.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.first.rangeJR(j).firstIndex(); k<=basis.first.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x2;
    for (int k= basis.second.mra.rangeI(basis.second.j0).firstIndex();
             k<=basis.second.mra.rangeI(basis.second.j0).lastIndex(); ++k) {
        LambdaB_x2.insert(Index1D(basis.second.j0,k,XBSpline));
    }
    for (int j=basis.second.j0; j<=J; ++j) {
        for (int k=basis.second.rangeJL(j).firstIndex(); k<=basis.second.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.second.rangeJR(j).firstIndex(); k<=basis.second.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x3;
    for (int k= basis.third.mra.rangeI(basis.third.j0).firstIndex();
             k<=basis.third.mra.rangeI(basis.third.j0).lastIndex(); ++k) {
        LambdaB_x3.insert(Index1D(basis.third.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.third.rangeJL(j).firstIndex(); k<=basis.third.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x3.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.third.rangeJR(j).firstIndex(); k<=basis.third.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x3.insert(Index1D(j,k,XWavelet));
        }
    }
    for (const_set1d_it it_x1=LambdaB_x1.begin(); it_x1!=LambdaB_x1.end(); ++it_x1) {
        for (const_set1d_it it_x2=LambdaB_x2.begin(); it_x2!=LambdaB_x2.end(); ++it_x2) {
            for (const_set1d_it it_x3=LambdaB_x3.begin(); it_x3!=LambdaB_x3.end(); ++it_x3) {
                Index3D newindex((*it_x1), (*it_x2), (*it_x3));
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
    }
}


template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end())  return;
    else                            v[index1d] = 0.;

    int  j = index1d.j;
    long k = index1d.k;

    Support<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                if (overlap(supp,new_supp)>0) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                }
            }
        }
        else {
            bool foundPredecessor = false;
            // First, we check whether there is a predecessor.
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)!=v.end()) {
                        foundPredecessor = true;
                        break;
                    }
                }
            }
            // Second, in case there exists no predecessor in the tree, we add the first candidate
            // to the tree.
            if (!foundPredecessor) {
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                        Index1D new_index1d(Index1D(new_j,new_k,new_type));
                        completeMultiTree(basis,new_index1d,v,sparsetree);
                        break;
                    }
                }
            }
        }
    }
    return;
}

template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v, int coordDirec, bool sparsetree,
                  bool isAlreadyMultiTree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (isAlreadyMultiTree) {
        if (v.find(index2d)!=v.end())  return;
        else                           v[index2d] = 0.;
    }
    else {
        //std::cerr << "     Completion to multitree for index = " << index2d << std::endl;
        if (v.find(index2d)==v.end())  v[index2d] = 0.;
        //std::cerr << "     Node inserted for index = " << index2d << std::endl;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
        Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
        //check x-direction
        int new_j_x = 0;
        long new_k_x_first = 0, new_k_x_last = 0;
        bool checkPredecessors=true;
        XType new_type_x = XWavelet;
        if (j_x==j0_x && index_x.xtype==XWavelet) {
            basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XBSpline;
            assert(new_j_x==j_x);
        }
        else if (j_x>j0_x && index_x.xtype==XWavelet) {
            basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XWavelet;
            assert(new_j_x==j_x-1);
        }
        else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (overlap(supp_x,new_supp_x)>0) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }
                if (!foundPredecessor) {
                    for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                        Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                        if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                            Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                            completeMultiTree(basis,new_index2d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec==0 || coordDirec==2) {
        Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
        //check y-direction
        bool checkPredecessors=true;
        int new_j_y = 0;
        long new_k_y_first = 0, new_k_y_last = 0;
        XType new_type_y = XWavelet;
        if (j_y==j0_y && index_y.xtype==XWavelet) {
            basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XBSpline;
            assert(new_j_y==j_y);
        }
        else if (j_y>j0_y && index_y.xtype==XWavelet) {
            basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XWavelet;
            assert(new_j_y==j_y-1);
        }
        else {
            checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (overlap(supp_y,new_supp_y)>0) {
                        Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                        Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }
                if (!foundPredecessor) {
                    for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                        Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                        if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                            Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
                            completeMultiTree(basis,new_index2d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v, int coordDirec, bool sparsetree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    if (v.find(index3d)!=v.end())  return;
    else                           v[index3d] = 0.;

    Index1D index_x = index3d.index1;
    Index1D index_y = index3d.index2;
    Index1D index_z = index3d.index3;

    int  j_x = index_x.j, j_y = index_y.j, j_z = index_z.j;
    long k_x = index_x.k, k_y = index_y.k, k_z = index_z.k;

    if (coordDirec==0 || coordDirec==1) {
        Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
        //check x-direction
        int new_j_x = 0;
        long new_k_x_first = 0, new_k_x_last = 0;
        bool checkPredecessors=true;
        XType new_type_x = XWavelet;
        if (j_x==j0_x && index_x.xtype==XWavelet) {
            basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XBSpline;
            assert(new_j_x==j_x);
        }
        else if (j_x>j0_x && index_x.xtype==XWavelet) {
            basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XWavelet;
            assert(new_j_x==j_x-1);
        }
        else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (overlap(supp_x,new_supp_x)>0) {
                        Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                        Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }

                if (!foundPredecessor) {
                    for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                        Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                        if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                            Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec==0 || coordDirec==2) {
        Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
        //check y-direction
        bool checkPredecessors=true;
        int new_j_y = 0;
        long new_k_y_first = 0, new_k_y_last = 0;
        XType new_type_y = XWavelet;
        if (j_y==j0_y && index_y.xtype==XWavelet) {
            basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XBSpline;
            assert(new_j_y==j_y);
        }
        else if (j_y>j0_y && index_y.xtype==XWavelet) {
            basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XWavelet;
            assert(new_j_y==j_y-1);
        }
        else {
            checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (overlap(supp_y,new_supp_y)>0) {
                        Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                        Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }

                if (!foundPredecessor) {
                    for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                        Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                        if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                            Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec== 0 || coordDirec==3) {
        Support<typename Basis::T> supp_z = basis.third.generator(index_z.xtype).support(j_z,k_z);
        //check z-direction
        bool checkPredecessors=true;
        int new_j_z = 0;
        long new_k_z_first = 0, new_k_z_last = 0;
        XType new_type_z = XWavelet;
        if (j_z==j0_z && index_z.xtype==XWavelet) {
            basis.third.getScalingNeighborsForWavelet(j_z,k_z,basis.third,new_j_z,new_k_z_first,new_k_z_last);
            new_type_z = XBSpline;
            assert(new_j_z==j_z);
        }
        else if (j_z>j0_z && index_z.xtype==XWavelet) {
            basis.third.getLowerWaveletNeighborsForWavelet(j_z,k_z,basis.third,new_j_z,new_k_z_first,new_k_z_last);
            new_type_z = XWavelet;
            assert(new_j_z==j_z-1);
        }
        else {
            checkPredecessors=false; // no "return" here, see above
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                    Support<typename Basis::T> new_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                    if (overlap(supp_z,new_supp_z)>0) {
                        Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                    Support<typename Basis::T> covered_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                    if (covered_supp_z.l1<=supp_z.l1 && covered_supp_z.l2>=supp_z.l2) {
                        Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }
                if (!foundPredecessor) {
                    for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                        Support<typename Basis::T> covered_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                        if (covered_supp_z.l1<=supp_z.l1 && covered_supp_z.l2>=supp_z.l2) {
                            Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2 && coordDirec !=3 ) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

/*
template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    if (v.find(index3d)!=v.end())  return;
    else                           v[index3d] = 0.;


    Index1D index_x = index3d.index1;
    Index1D index_y = index3d.index2;
    Index1D index_z = index3d.index3;

    int  j_x = index_x.j, j_y = index_y.j, j_z = index_z.j;
    long k_x = index_x.k, k_y = index_y.k, k_z = index_z.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
    //check x-direction
    if (j_x==j0_x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(j_x==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index3D new_index3d(Index1D(j_scaling,k_scaling,XBSpline),index_y,index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_x>j0_x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_x-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index3D new_index3d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y,index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
    //check y-direction
    if (j_y==j0_y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(j_y==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index3D new_index3d(index_x,Index1D(j_scaling,k_scaling,XBSpline),index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_y>j0_y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_y-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index3D new_index3d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet),index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }

    Support<typename Basis::T> supp_z = basis.third.generator(index_z.xtype).support(j_z,k_z);
    //check z-direction
    if (j_z==j0_z  && index_z.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.third.getScalingNeighborsForWavelet(j_z,k_z,basis.third,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(j_z==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_z = basis.third.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_z,new_supp_z)>0) {
                Index3D new_index3d(index_x,index_y,Index1D(j_scaling,k_scaling,XBSpline));
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_z>j0_z && index_z.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.third.getLowerWaveletNeighborsForWavelet(j_z,k_z,basis.third,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_z-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_z = basis.third.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_z,new_supp_z)>0) {
                Index3D new_index3d(index_x,index_y,Index1D(j_wavelet,k_wavelet,XWavelet));
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    return;
}
*/




template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index2D> &v, int j, T gamma)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    for (long k1=basis.first.mra.rangeI(j0_x).firstIndex(); k1<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
        for (long k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
            Index1D row(j0_x,k1,XBSpline);
            Index1D col(j0_y,k2,XBSpline);
            v[Index2D(row,col)] = 0.;
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0_y+i2-1;
            for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_x,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                v[Index2D(row,col)] = 0.;
                v[Index2D(col,row)] = 0.;
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_x+i1-1;
        for (int i2=1; i2<=j; ++i2) {
            if (T(i1+i2)-gamma*std::max(i1,i2)>(1-gamma)*j) {
                continue;
            }
            int j2=j0_y+i2-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    v[Index2D(row,col)] = 0.;
                }
            }
        }
    }
    return;
}

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index3D> &v, int J, T gamma)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    for (int i_x=0; i_x<=J; ++i_x) {
        Range<int> range_x(0,0);
        XType xtype_x;
        int j_x;
        if (i_x==0) {   range_x = basis.first.mra.rangeI(j0_x);   j_x = j0_x;       xtype_x = XBSpline; }
        else        {   range_x = basis.first.rangeJ(j0_x+i_x-1); j_x = j0_x+i_x-1; xtype_x = XWavelet; }
        for (int k_x=range_x.firstIndex(); k_x<=range_x.lastIndex(); ++k_x) {
            Index1D index_x(j_x, k_x, xtype_x);

            for (int i_y=0; i_y<=J-i_x; ++i_y) {
                Range<int> range_y(0,0);
                XType xtype_y;
                int j_y;
                if (i_y==0) {   range_y = basis.second.mra.rangeI(j0_y);   j_y = j0_y;       xtype_y = XBSpline; }
                else        {   range_y = basis.second.rangeJ(j0_y+i_y-1); j_y = j0_y+i_y-1; xtype_y = XWavelet; }
                for (int k_y=range_y.firstIndex(); k_y<=range_y.lastIndex(); ++k_y) {
                    Index1D index_y(j_y, k_y, xtype_y);

                    for (int i_z=0; i_z<=J-i_x-i_y; ++i_z) {
                        int sum = i_x+i_y+i_z;
                        int max = std::max(i_x,std::max(i_y,i_z));
                        if (sum - gamma*max > (1-gamma)*J) continue;
                        Range<int> range_z(0,0);
                        XType xtype_z;
                        int j_z;
                        if (i_z==0) {   range_z = basis.second.mra.rangeI(j0_z);   j_z = j0_z;       xtype_z = XBSpline; }
                        else        {   range_z = basis.second.rangeJ(j0_z+i_z-1); j_z = j0_z+i_z-1; xtype_z = XWavelet; }
                        for (int k_z=range_z.firstIndex(); k_z<=range_z.lastIndex(); ++k_z) {
                            Index1D index_z(j_z, k_z, xtype_z);
                            v[Index3D(index_x,index_y,index_z)] = 0.;
                        }
                    }
                }
            }
        }
    }
    return;
}


template <typename Basis>
void
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (jx==j0x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(jx==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jx>j0x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jx-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (jy==j0y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(jy==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jy>j0y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jy-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
}

template <typename Basis>
void
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (index_x.j==j0_x) {
        for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.first.rangeJ(index_x.j-1).firstIndex(),(long)index_x.k/2 - offset);
        long k_last  = std::min((long)basis.first.rangeJ(index_x.j-1).lastIndex(),(long)index_x.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (index_y.j==j0_y) {
        for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.second.rangeJ(index_y.j-1).firstIndex(),(long)index_y.k/2 - offset);
        long k_last  = std::min((long)basis.second.rangeJ(index_y.j-1).lastIndex(),(long)index_y.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
}

}   // namespace lawa

/*
 *
bool
comp_support(const Support<double> &supp1, const Support<double> &supp2) {
return supp1.l1<=supp2.l1 ? true : false;
}
 *
 *
std::vector<Support<T> > intersectingSupports;
bool foundPredecessor = false;
for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
    Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
    if (v.find(new_index3d)!=v.end()) {
        Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
        intersectingSupports.push_back(covered_supp_x);
        if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
            foundPredecessor = true;
            break;
        }
    }
}

if (!foundPredecessor && intersectingSupports.size()>0) {
    std::cerr << "Before sorting: " << std::endl;
    for (typename std::vector<Support<T> >::const_iterator it=intersectingSupports.begin();
          it!=intersectingSupports.end(); ++it) {
        std::cerr << (*it) << " " << supp_x <<  std::endl;
    }
    std::cerr << "After sorting:  " << std::endl;
    sort(intersectingSupports.begin(), intersectingSupports.end(), comp_support);
    for (typename std::vector<Support<T> >::const_iterator it=intersectingSupports.begin();
          it!=intersectingSupports.end(); ++it) {
        std::cerr << (*it) << " " << supp_x << std::endl;
    }
    std::cerr << std::endl;
    typename std::vector<Support<T> >::const_iterator it1=intersectingSupports.begin();
    typename std::vector<Support<T> >::const_iterator it2=intersectingSupports.begin();
    ++it2;
    Support<T> jointSupport = (*it1);
    while (it2!=intersectingSupports.end()) {
        if (overlap((*it1),(*it2))) {
            jointSupport.l1 = std::min((*it1).l1, jointSupport.l1);
            jointSupport.l2 = std::max((*it2).l2, jointSupport.l2);
        }
        else {
            break;
        }
        ++it1; ++it2;
    }
    if (jointSupport.l1<=supp_x.l1 && jointSupport.l2>=supp_x.l2) foundPredecessor=true;
}
 */
