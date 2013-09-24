//
//  adaptiveoperator2d.h
//  LAWA
//
//  Created by Kristina Steih on 29.03.12.
//  Copyright 2012 Universität Ulm. All rights reserved.
//

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>

namespace lawa {
  
  template <typename T>
  struct AdaptiveOperator2D : Operator2D<T> {
    
    typedef flens::SparseGeMatrix<CRS<T,CRS_General> >                  SparseMatrixT;

    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) = 0;
    
    virtual Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x) = 0;
    
    virtual void
    toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol, 
												SparseMatrixT &A, T eps, bool useLinearIndex=false) = 0; 
    
  };
  
} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H
