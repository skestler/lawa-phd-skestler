/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


namespace lawa {


template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens)
{
    //std::cout << "   Assembling of sparse matrix started..." << std::endl;
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    std::map<Index,int,lt<Lexicographical,Index> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    A.compression.setParameters(LambdaRow);
    for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        IndexSet<Index> LambdaRowSparse = A.compression.SparsityPattern(*col, LambdaRow);
        for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
            T tmp = A(*row,*col);
            if (fabs(tmp)>0)                A_flens(row_indices[*row],col_count) = tmp;
        }
    }
    A_flens.finalize();
    //std::cout << "   N^2 = " << LambdaRow.size()*LambdaCol.size() << ", nonzeros = " << A_flens.numNonZeros() << std::endl;
}

template <typename T, typename Index, typename SpaceIndex, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRowOp, const IndexSet<SpaceIndex>& LambdaRowInitCond,
                    const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens)
{
    std::cout << "   Assembling of sparse matrix started..." << std::endl;
    typedef typename IndexSet<Index>::const_iterator const_set_op_it;
    typedef typename IndexSet<SpaceIndex>::const_iterator const_set_initcond_it;
    std::map<Index,int,lt<Lexicographical,Index> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_op_it row=LambdaRowOp.begin(); row!=LambdaRowOp.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    A.c.setParameters(LambdaRowOp);

    for (const_set_op_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        IndexSet<Index> LambdaRowSparse = A.c.SparsityPattern(*col, LambdaRowOp);
        for (const_set_op_it row=LambdaRowOp.begin(); row!=LambdaRowOp.end(); ++row) {
            T tmp = A(*row,*col);
            if (fabs(tmp)>0)                A_flens(row_indices[*row],col_count) = tmp;
        }
    }
    col_count = 1;

    for (const_set_op_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        row_count = LambdaRowOp.size()+1;
        for (const_set_initcond_it row=LambdaRowInitCond.begin(); row!=LambdaRowInitCond.end(); ++row, ++row_count) {
            T tmp = A(*row, *col);
            if (fabs(tmp)>0)                A_flens(row_count,col_count) = tmp;
        }
    }
    A_flens.finalize();
    std::cout << "   ...finished." << std::endl;
}

template <typename T, typename RowIndex, typename ColIndex, typename MA>
Coefficients<Lexicographical,T,RowIndex>
mv(const IndexSet<RowIndex> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,ColIndex > &v)
{
    typedef typename IndexSet<RowIndex>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,ColIndex>::const_iterator coeff_const_it;

    Coefficients<Lexicographical,T,RowIndex> w;
    Timer timer;
    timer.start();
    for (set_const_it lambda = LambdaRow.begin(); lambda != LambdaRow.end(); ++lambda) {
        T val = 0;
        for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
            val += A(*lambda,(*mu).first) * (*mu).second;
        }
        w[*lambda] = val;
    }
    timer.stop();
    //std::cout << "   Elapsed time for standard mv: " << timer.elapsed() << std::endl;

    return w;
}

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv_sparse(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    typedef typename IndexSet<Index>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

    Timer timer;
    timer.start();
    A.compression.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index> w_;
    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
        IndexSet<Index> LambdaRowSparse = A.compression.SparsityPattern((*mu).first, LambdaRow);
        for (set_const_it lambda=LambdaRowSparse.begin(); lambda!=LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += A(*lambda,(*mu).first) * (*mu).second;
        }
    }
    timer.stop();
    //std::cout << "   Elapsed time for improved mv: " << timer.elapsed() << std::endl;

    return w_;
}

/*
template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    // ersetzen durch iteration Ÿber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
    typedef typename IndexSet<Index>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

    Coefficients<Lexicographical,T,Index> w(v.d,v.d_);
    for (set_const_it lambda = LambdaRow.begin(); lambda != LambdaRow.end(); ++lambda) {
        T val = 0;
        for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
            val += A(t, *lambda,(*mu).first) * (*mu).second;
        }
        w[*lambda] = val;
    }
    return w;
}

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv_sparse(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator coeff_const_it;

    A.c.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index2D> w_(v.d,v.d_);
    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
        IndexSet<Index2D> LambdaRowSparse = A.c.SparsityPattern((*mu).first, LambdaRow);
        for (set2d_const_it lambda = LambdaRowSparse.begin(); lambda != LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += A(t,*lambda,(*mu).first) * (*mu).second;
        }
    }
    return w_;
}
*/

/*
template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
         const Coefficients<Lexicographical,T,Index > &f, T &residual, T tol,
         bool useOptimizedAssembling, int maxIterations)
{
    typedef typename IndexSet<Index >::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

    Timer timer;
    int N = Lambda.size();
    flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
    if (useOptimizedAssembling) {
        A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
    }
    else {
        toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
    }

    if (Lambda.size() > 0) {
      std::cout << "    Build Dense Vectors ..." << std::endl;
      timer.start();
        DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
        int row_count=1;
        const_coeff_it f_end = f.end();
        const_coeff_it u_end = u.end();
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
            const_coeff_it f_it = f.find(*row);
            if (f_it != f_end) rhs(row_count) = (*f_it).second;
            else               rhs(row_count) = 0.;
            const_coeff_it u_it = u.find(*row);
            if (u_it != u_end) x(row_count) = (*u_it).second;
            else               x(row_count) = 0.;

        }
      timer.stop();
      std::cout << "    .... done : " << timer.elapsed() << " seconds " << std::endl;  

      std::cout << "    Start cg ... " << std::endl;
      timer.start();
        int number_of_iterations = lawa::cg(A_flens,x,rhs, tol, maxIterations);
      timer.stop();
            
      std::cout << "    .... done : " << timer.elapsed() << " seconds " << std::endl;
        Ax = A_flens*x;
        res= Ax-rhs;
        residual = std::sqrt(res*res);
        row_count = 1;
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
            //const_coeff_it u_it = u.find(*row);
            //if (u_it != u_end) u[*row] = x(row_count);
            //else               u[*row] = 0.;
            u[*row] = x(row_count);
        }
        
        return number_of_iterations;
    }
    else return -1;

}

template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
            const Coefficients<Lexicographical,T,Index > &f, T &residual, T tol,
            bool useOptimizedAssembling, int maxIterations)
{
        typedef typename IndexSet<Index >::const_iterator const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

        int N = Lambda.size();
        flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
        if (useOptimizedAssembling) {
            A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
        }
        else {
            toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
        }

        if (Lambda.size() > 0) {
            DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
            int row_count=1;
            const_coeff_it f_end = f.end();
            const_coeff_it u_end = u.end();
            for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
                const_coeff_it f_it = f.find(*row);
                if (f_it != f_end) rhs(row_count) = (*f_it).second;
                else               rhs(row_count) = 0.;
                const_coeff_it u_it = u.find(*row);
                if (u_it != u_end) x(row_count) = (*u_it).second;
                else               x(row_count) = 0.;

            }
            int number_of_iterations = lawa::gmres(A_flens,x,rhs, tol, maxIterations);
            Ax = A_flens*x;
            res= Ax-rhs;
            residual = std::sqrt(res*res);
            row_count = 1;
            for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
                u[*row] = x(row_count);
            }
            return number_of_iterations;
        }
        else return -1;

}

template <typename T, typename Index, typename MA>
int
GMRESM_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
            const Coefficients<Lexicographical,T,Index > &f, T &residual, T tol,
            bool useOptimizedAssembling, int maxIterations, int m)
{
    typedef typename IndexSet<Index >::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;
    
    int N = Lambda.size();
    flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
    if (useOptimizedAssembling) {
        A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
    }
    else {
        toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
    }
    
    if (Lambda.size() > 0) {
        DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
        int row_count=1;
        const_coeff_it f_end = f.end();
        const_coeff_it u_end = u.end();
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
            const_coeff_it f_it = f.find(*row);
            if (f_it != f_end) rhs(row_count) = (*f_it).second;
            else               rhs(row_count) = 0.;
            const_coeff_it u_it = u.find(*row);
            if (u_it != u_end) x(row_count) = (*u_it).second;
            else               x(row_count) = 0.;
            
        }
        int number_of_iterations = lawa::gmresm(A_flens,x,rhs, tol, m, maxIterations);
        Ax = A_flens*x;
        res= Ax-rhs;
        residual = std::sqrt(res*res);
        row_count = 1;
        for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
            u[*row] = x(row_count);
        }
        return number_of_iterations;
    }
    else return -1;
    
}
*/

/*
template <typename T, typename Index, typename MA>
int
GMRES_Solve_PG(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol,
           MA &A, Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f, T &residual, T tol, int maxIterations)
{
		if(LambdaRow.size() != LambdaCol.size()){
			std::cerr << "Cannot apply GMRES: Matrix not quadratic! " << std::endl;
			exit(1);
		}
		
        typedef typename IndexSet<Index >::const_iterator const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

        int NumOfRows = (int)LambdaRow.size();
        int NumOfCols = (int)LambdaCol.size();
        flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
        toFlensSparseMatrix(A, LambdaRow, LambdaCol, A_flens);
    
  
    //----------------------------
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    densify(cxxblas::NoTrans, A_flens, A_dense);
    std::ofstream A_matrix("LHS_Matrix.txt");
    A_matrix << A_dense << std::cout;
    A_matrix.close();
    //----------------------------


        if (LambdaRow.size() > 0) {
            DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), res(NumOfRows), Ax(NumOfRows);
            int row_count=1;
            for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
                if (f.count((*row)) > 0) {
                    const_coeff_it it = f.find(*row);
                    rhs(row_count) = (*it).second;
                }
                else                     rhs(row_count) = 0.;
            }
            
            //----------------------------
            std::ofstream b_vector("RHS_Vector.txt");
            b_vector << rhs << std::cout;
            b_vector.close();
            //----------------------------
            
            std::cout << "Starting gmres (pg)..." << std::endl;
            int number_of_iterations = lawa::gmres(A_flens,x,rhs, tol, maxIterations);
            std::cout << "... finished" << std::endl;
            Ax = A_flens*x;
            res= Ax-rhs;
            residual = std::sqrt(res*res);
            int col_count = 1;
            for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
                u[*col] = x(col_count);
            }
            return number_of_iterations;
        }
        else return -1;

}

template <typename T, typename Index, typename SpaceIndex, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRowOp, const IndexSet<SpaceIndex> &LambdaRowInitCond,
           MA &A, const IndexSet<Index> &LambdaCol,
           Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f,
           const Coefficients<Lexicographical,T,SpaceIndex > &u0,
           T &r, T tol, int maxIterations)
{
    typedef typename IndexSet<Index >::const_iterator const_set_op_it;
    typedef typename IndexSet<SpaceIndex >::const_iterator const_set_initcond_it;
    typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,SpaceIndex >::const_iterator const_coeff_initcond_it;
    typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

    std::cerr << "CGLS_SOLVE called..." << std::endl;
    int NumOfCols = LambdaCol.size();
    int NumOfRows = LambdaRowOp.size() + LambdaRowInitCond.size();
    flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
    toFlensSparseMatrix(A, LambdaRowOp, LambdaRowInitCond, LambdaCol, A_flens);

    if (LambdaCol.size() > 0) {
        DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), res(NumOfRows), Ax(NumOfRows);
        int row_count=1;
        const_coeff_it f_end = f.end();
        for (const_set_op_it row=LambdaRowOp.begin(); row!=LambdaRowOp.end(); ++row, ++row_count) {
            const_coeff_it f_it = f.find(*row);
            if (f_it != f_end) rhs(row_count) = (*f_it).second;
            else               rhs(row_count) = 0.;
        }

        for (const_set_initcond_it row=LambdaRowInitCond.begin(); row!=LambdaRowInitCond.end(); ++row, ++row_count) {
            if (u0.count((*row)) > 0) {
                const_coeff_initcond_it it = u0.find(*row);
                rhs(row_count) = (*it).second;
            }
            else                      rhs(row_count) = 0.;
        }
        const_coeff_it u_end = u.end();
        int col_count=1;
        for (const_set_op_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
            const_coeff_it u_it = u.find(*col);
            if (u_it != u_end) x(col_count) = (*u_it).second;
            else               x(col_count) = 0.;
        }
        int number_of_iterations = lawa::cgls(A_flens,x,rhs, tol, maxIterations);
        Ax = A_flens*x;
        res= Ax-rhs;
        r = std::sqrt(res*res);
        row_count = 1;
        for (const_set_op_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row, ++row_count) {
            u[*row] = x(row_count);
        }
        return number_of_iterations;
    }
    else return -1;
}

template <typename T, typename Index, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol,
           MA &A, Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f, T &r, T tol, int maxIterations)
{
        typedef typename IndexSet<Index >::const_iterator const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

        int NumOfRows = (int)LambdaRow.size();
        int NumOfCols = (int)LambdaCol.size();
        flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
        toFlensSparseMatrix(A, LambdaRow, LambdaCol, A_flens);

        if (LambdaRow.size() > 0) {
            DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), res(NumOfRows), Ax(NumOfRows);
            int row_count=1;
            for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
                if (f.count((*row)) > 0) {
                    const_coeff_it it = f.find(*row);
                    rhs(row_count) = (*it).second;
                }
                else                     rhs(row_count) = 0.;
            }
            std::cout << "Starting cgls..." << std::endl;
            int number_of_iterations = lawa::cgls(A_flens,x,rhs, tol, maxIterations);
            std::cout << "... finished" << std::endl;
            Ax = A_flens*x;
            res= Ax-rhs;
            r= std::sqrt(res*res);
            int col_count = 1;
            for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
                u[*col] = x(col_count);
            }
            return number_of_iterations;
        }
        else return -1;

}
*/
}   //namespace lawa

/*
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    densify(cxxblas::NoTrans,A_flens,A_dense);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(NumOfRows,NumOfRows),V(NumOfCols,NumOfCols);
    DenseVector<Array<T> > s(NumOfCols);
    std::cout << "Computing svd..." << std::endl;
    int iterations = svd(A_dense,s,U,V);
    std::cout << " ... finished after " << iterations << std::endl;
    std::cout << "Largest singular value: " << s(s.firstIndex()) << ", smallest singular value: " << s(s.lastIndex()) << std::endl;
*/

/*
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense1;
    densify(cxxblas::NoTrans,A_flens,A_dense1);
    //std::cout << A_dense1 << std::endl;
    DenseVector<Array<T> > wr(N), wi(N);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > vl,vr;
    ev(false, false, A_dense1, wr, wi, vl, vr);
    T cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
    std::cout << "Largest eigenvalue: " << CB << ", smallest eigenvalue: " << cB << std::endl;
*/

