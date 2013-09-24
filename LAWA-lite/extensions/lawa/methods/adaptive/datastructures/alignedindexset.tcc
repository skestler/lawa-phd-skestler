namespace lawa {

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::AlignedIndexSet(void)
: map(), n1(0), n2(0)
{

}

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::AlignedIndexSet(size_t _n1, size_t _n2)
: map(_n1), n1(_n1), n2(_n2)
{

}

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::align_x1(const IndexSet<Index> &Lambda) {

    for (const_set_index_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        //todo: works only in 2d!! additional routine for index splitting required.
        map_prinindex_it p_prinindex=map.find((*it).index1);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.insert((*it).index2);
        }
        else {
            IndexSet<PrincipalIndex> Lambda_x2(n2);
            Lambda_x2.insert((*it).index2);
            map[(*it).index1] = Lambda_x2;
        }
    }
    return;
}

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::align_x2(const IndexSet<Index> &Lambda) {

    for (const_set_index_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        //todo: works only in 2d!! additional routine for index splitting required.
        map_prinindex_it p_prinindex=map.find((*it).index2);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.insert((*it).index1);
        }
        else {
            IndexSet<PrincipalIndex> Lambda_x1(n2);
            Lambda_x1.insert((*it).index1);
            map[(*it).index2] = Lambda_x1;
        }
    }
    return;
}

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::clear()
{
    for (map_prinindex_it it=map.begin(); it!=map.end(); ++it) {
        (*it).second.clear();
    }
    map.clear();
}

template <typename Index, typename PrincipalIndex, typename AlignedIndex>
std::ostream& operator<< (std::ostream &s,
                          const AlignedIndexSet<Index,PrincipalIndex,AlignedIndex> &alignedLambda)
{
    s << std::endl << "AlignedIndexSet:" << std::endl;
    for (typename AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::const_map_prinindex_it row=alignedLambda.map.begin();
            row!=alignedLambda.map.end(); ++row)
    {
        for (typename AlignedIndexSet<Index,PrincipalIndex,AlignedIndex>::const_set_aligindex_it col=(*row).second.begin();
                col!=(*row).second.end(); ++col)
        {
            s << (*row).first << ", " << *col << std::endl;
        }
    }
    return s << std::endl;
}

}   // namespace lawa
