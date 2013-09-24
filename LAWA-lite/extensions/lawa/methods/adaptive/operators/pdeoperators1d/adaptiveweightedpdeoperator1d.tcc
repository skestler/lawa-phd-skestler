namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
AdaptiveWeightedPDEOperator1D<T,Side,Domain,Cons>::
AdaptiveWeightedPDEOperator1D(const Basis1D& _basis1d, Function<T> &_reaction_f,
                              Function<T> &_convection_f, Function<T>& _diffusion_f, int order,
                              bool reactionIsZero, bool convectionIsZero,
                              bool diffusionIsZero)
: basis1d(_basis1d),
  compression1d(basis1d),
  weightedpdeop1d(basis1d,_reaction_f,_convection_f,_diffusion_f,order,reactionIsZero,convectionIsZero,diffusionIsZero),
  prec1d(), data(weightedpdeop1d, prec1d, compression1d)
{

}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveWeightedPDEOperator1D<T,Side,Domain,Cons>::operator()(const Index1D &row_index,
                                                              const Index1D &col_index)
{
    return data(row_index, col_index);
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
AdaptiveWeightedPDEOperator1D<T,Side,Domain,Cons>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    return this->operator()(row_index,col_index);
}

}   //namespace lawa
