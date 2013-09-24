namespace lawa {

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection CoordX>
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex,CoordX>::
AlignedCoefficients(void)
: map(), n1(0), n2(0)
{

}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection CoordX>
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex,CoordX>::
AlignedCoefficients(size_t _n1,size_t _n2)
: map(_n1), n1(_n1), n2(_n2)
{

}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection CoordX>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex,CoordX>
::align(const Coefficients<Lexicographical,T,Index> &coeff, short J)
{
    Split<Index,PrincipalIndex,AlignedIndex,CoordX> split;
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        PrincipalIndex prinIndex;
        AlignedIndex   aligIndex;
        split((*it).first,prinIndex,aligIndex);
        map_prinindex_it p_prinindex=map.find(prinIndex);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
        else {
            size_t tmp = std::max((size_t)pow2i<long>(J-prinIndex.levelSum()+2),n2);
            Coefficients<Lexicographical,T,AlignedIndex> coeff_x2;
            //coeff_x2[aligIndex] = (*it).second;
            map[prinIndex] = coeff_x2;
            map_prinindex_it p_prinindex=map.find(prinIndex);
            (*p_prinindex).second.Rehash(tmp);
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
    }
    //std::cerr << "AlignedCoefficients: " << map.size() << " " << map.bucket_count() << std::endl;
}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection CoordX>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex,CoordX>
::align(const Coefficients<Lexicographical,T,Index> &coeff, IndexSet<PrincipalIndex> &prinIndices,
        short J)
{
    Split<Index,PrincipalIndex,AlignedIndex,CoordX> split;
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        PrincipalIndex prinIndex;
        AlignedIndex   aligIndex;
        split((*it).first,prinIndex,aligIndex);
        map_prinindex_it p_prinindex=map.find(prinIndex);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
        else {
            prinIndices.insert(prinIndex);
            size_t tmp = std::max((size_t)pow2i<long>(J-prinIndex.levelSum()+2),n2);
            Coefficients<Lexicographical,T,AlignedIndex> coeff_x2;
            //coeff_x2[aligIndex] = (*it).second;
            map[prinIndex] = coeff_x2;
            map_prinindex_it p_prinindex=map.find(prinIndex);
            (*p_prinindex).second.Rehash(tmp);
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
    }
}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex,
          CoordinateDirection CoordX>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex,CoordX>
::align_ExcludeAndOrthogonal(const Coefficients<Lexicographical,T,Index> &coeff,
                             const Coefficients<Lexicographical,T,Index> &exclude,
                             const IndexSet<PrincipalIndex> &orthogonal, short J)
{
    Split<Index,PrincipalIndex,AlignedIndex,CoordX> split;
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        PrincipalIndex prinIndex;
        AlignedIndex   aligIndex;
        split((*it).first,prinIndex,aligIndex);
        if (orthogonal.find(prinIndex)==orthogonal.end()) continue;

        if (exclude.find((*it).first)!=exclude.end()) continue;

        map_prinindex_it p_prinindex=map.find(prinIndex);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
        else {
            size_t tmp = std::max((size_t)pow2i<long>(J-prinIndex.levelSum()+2),n2);
            Coefficients<Lexicographical,T,AlignedIndex> coeff_x2;
            //coeff_x2[aligIndex] = (*it).second;
            map[prinIndex] = coeff_x2;
            map_prinindex_it p_prinindex=map.find(prinIndex);
            (*p_prinindex).second.Rehash(tmp);
            (*p_prinindex).second.operator[](aligIndex) = (*it).second;
        }
    }
}

/*
template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
AlignedCoefficients2<T,Index,PrincipalIndex,AlignedIndex>::AlignedCoefficients2(void)
: principalIndices(24593), principalIndexToAlignedIndices()
{

}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedCoefficients2<T,Index,PrincipalIndex,AlignedIndex>
::align_x1(const Coefficients<Lexicographical,T,Index> &coeff)
{
    int n=0;
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        const_coeff_prinindex_it prin_it = principalIndices.find((*it).first.index1);
        if (prin_it!=principalIndices.end()) {
            int pos = (*prin_it).second;
            const AlignedIndex* tmp = &((*it).first.index2);
            principalIndexToAlignedIndices[pos].push_back(tmp);
        }
        else {
            principalIndices[(*it).first.index1] = n;
            AlignedIndices alignedIndices;
            const AlignedIndex* tmp = &((*it).first.index2);
            alignedIndices.push_back(tmp);
            principalIndexToAlignedIndices.push_back(alignedIndices);
//            principalIndexToAlignedIndices[n].push_back(tmp);
            ++n;
        }
    }
    //std::cerr << "AlignedCoefficients: " << map.size() << " " << map.bucket_count() << std::endl;
}
*/

}   // namespace lawa
