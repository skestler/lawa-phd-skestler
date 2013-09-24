namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis, QuadratureType Quad>
WeightedConvectionOperator1D_PG<T,TrialBasis,TestBasis,Quad>::
WeightedConvectionOperator1D_PG(const TrialBasis& _trialbasis,
		const TestBasis& _testbasis, Function<T> weightFct)
    : trialbasis(_trialbasis), testbasis(_testbasis), W(weightFct), integral(W, _testbasis, _trialbasis)
{
}

template <typename T, typename TrialBasis, typename TestBasis, QuadratureType Quad>
T
WeightedConvectionOperator1D_PG<T,TrialBasis,TestBasis,Quad>::
operator()(XType xtype1, int j1, int k1,
		   XType xtype2, int j2, int k2) const
{   
    // v * u_x
    return integral(j1, k1, xtype1, 0, j2, k2, xtype2, 1);
}

template <typename T, typename TrialBasis, typename TestBasis, QuadratureType Quad>
T
WeightedConvectionOperator1D_PG<T,TrialBasis,TestBasis,Quad>::
operator()(const Index1D &row_index,
		   const Index1D &col_index) const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                            col_index.xtype, col_index.j, col_index.k);
}

}    //namespace lawa
