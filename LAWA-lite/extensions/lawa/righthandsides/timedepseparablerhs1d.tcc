namespace lawa {
    
template<typename T, typename Basis1D>
TimedepSeparableRHS1D<T, Basis1D>::TimedepSeparableRHS1D(const Basis1D& _basis,
                                                         const SeparableFunction2D<T>& _F,
                                                         int order)
    : basis(_basis), F(_F), integralf(_F.F_y, basis)
{
    integralf.quadrature.setOrder(order);
}
    
template<typename T, typename Basis1D>
T
TimedepSeparableRHS1D<T, Basis1D>::operator()(T t, XType xtype, int j, long k) const
{
    return F.F_x(t) * integralf(j,k,xtype,0);
}   
   

template<typename T, typename Basis1D>
T
TimedepSeparableRHS1D<T, Basis1D>::operator()(T t, const Index1D &index) const
{
    return this->operator()(t, index.xtype, index.j, index.k);
} 

//==========================================================================================//

template<typename T, typename TimedepRHS>
SumOfTimedepRHS1D<T, TimedepRHS>::SumOfTimedepRHS1D(const TimedepRHS &_rhs1,
                                                    const TimedepRHS &_rhs2)
    : rhs1(_rhs1), rhs2(_rhs2)
{
}


template<typename T, typename TimedepRHS>
T
SumOfTimedepRHS1D<T, TimedepRHS>::operator()(T t, XType xtype, int j, long k) const
{
       return rhs1(t, xtype, j, k)
          +rhs2(t, xtype, j, k); 
}

template<typename T, typename TimedepRHS>
T
SumOfTimedepRHS1D<T, TimedepRHS>::operator()(T t, const Index1D &index) const
{
    return rhs1(t, index.xtype, index.j, index.k)
          +rhs2(t, index.xtype, index.j, index.k);
}

template<typename T, typename TimedepRHS>
Coefficients<Lexicographical,T,Index1D>
SumOfTimedepRHS1D<T, TimedepRHS>::operator()(T t, const IndexSet<Index1D> &Lambda) const
{
    typedef typename IndexSet<Index1D>::iterator const_set_it;
    Coefficients<Lexicographical,T,Index1D> ret;
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = this->operator()(t, *lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}

} // namespace lawa
