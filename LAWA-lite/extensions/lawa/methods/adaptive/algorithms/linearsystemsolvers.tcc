namespace lawa {

#define CGMRES_ITER 10

template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
         const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
         T &timeMatrixVector, int assemble_matrix)
{
    typedef typename IndexSet<Index >::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

    Timer timer;
    if (assemble_matrix==0) {

        T alpha, beta, rNormSquare, rNormSquarePrev;
        //Coefficients<Lexicographical,T,Index> r(SIZEHASHINDEX2D), p(SIZEHASHINDEX2D),
        //                                      Ap(SIZEHASHINDEX2D);
        size_t hms = u.size();
        Coefficients<Lexicographical,T,Index> r(hms), p(hms), Ap(hms);
        A.apply(u, tol/3.,Lambda, r);
        r -= f;
        p -= r;
        rNormSquare = r*r;

        for (int k=1; k<=maxIterations; k++) {
            if (sqrt(rNormSquare)<=tol) {
                res = sqrt(rNormSquare);
                return k;
            }
            Ap.setToZero();
            timer.start();
            A.apply(p,tol/3.,Lambda, Ap);
            timer.stop();
            timeMatrixVector = timer.elapsed();
            T pAp = p * Ap;
            alpha = rNormSquare/pAp;
            u += alpha*p;
            r += alpha*Ap;

            rNormSquarePrev = rNormSquare;
            rNormSquare = r*r;
            beta = rNormSquare/rNormSquarePrev;
            p *= beta;
            p -= r;
        }
        return maxIterations;
    }
    else {

        int N = Lambda.size();
        flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
        if (assemble_matrix==2) {
            A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
        }
        else {
            toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
        }

        if (Lambda.size() > 0) {
          std::cout << "    Build Dense Vectors ..." << std::endl;
          timer.start();
            DenseVector<Array<T> > rhs(N), x(N), residual(N), Ax(N);
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
            residual = Ax-rhs;
            res = std::sqrt(residual*residual);
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
    return -1;
}

template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
            int assemble_matrix)
{
        typedef typename IndexSet<Index >::const_iterator const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;
        typedef DenseVector<Array<T> >              DeVector;
        typedef GeMatrix<FullStorage<T, cxxblas::ColMajor> > DeMatrix;

        std::cerr << "GMRES_Solve called." << std::endl;

        int N = Lambda.size();
        if (assemble_matrix==0) {   //Algorithm by Saad, p.159

            Coefficients<Lexicographical,T,Index> v(2*u.size()), w_j(2*u.size());
            Coefficients<Lexicographical,T,Index> V[CGMRES_ITER+1];

            /*
            Coefficients<Lexicographical,T,Index> v(2*u.size()), w_j(2*u.size());
            Coefficients<Lexicographical,T,Index> V[CGMRES_ITER+1];
            for (int i=0;i<=CGMRES_ITER; ++i) {
                V[i].rehash(2*u.size());
            }
            */

            for (int iter=1; iter<=maxIterations/CGMRES_ITER; ++iter) {
                v.setToZero();
                for (int i=0; i<CGMRES_ITER+1; ++i) {
                    V[i].setToZero();
                }

                DeMatrix H(CGMRES_ITER+1,CGMRES_ITER);
                DeVector g(CGMRES_ITER+1);
                DeVector c(CGMRES_ITER+1), s(CGMRES_ITER+1);
                T nu, rho = tol + 1;
                T Htemp, h_ij, beta;


                A.apply(u, tol/3.,Lambda, v);
                v -= f;
                beta = v.norm(2.);
                std::cerr << "      GMRES_solve: current error = " << beta << std::endl;
                if (beta<tol) {
                    std::cerr << "      GMRES_solve: abort criterion met" << std::endl;
                    res = beta;
                    return (iter-1)*CGMRES_ITER;
                }

                v *= (-1./beta);
                V[0] = v;
                g(1) = beta;

                int j;
                for (j=0; j<CGMRES_ITER;) {
                    if (rho <= tol) {
                        break;
                    } else {
                           ++j;
                    }
                    w_j.setToZero();
                    A.apply(V[j-1], tol/3.,Lambda, w_j);
                    T normInitialWj = w_j.norm(2.);

                    for (int i=1; i<=j; ++i) {
                        H(i,j) = w_j * V[i-1];
                        w_j -= H(i,j) * V[i-1];
                    }
                    H(j+1,j) = w_j.norm(2.);
                    /*
                    if (H(j+1,j) / normInitialWj < 1.0) {
                        for (int i=1; i<=j; ++i) {
                            Htemp = w_j * V[i-1];
                            w_j -= Htemp * V[i-1];
                        }
                        H(j+1, j) = w_j.norm(2.);
                    }
                    */
                    for (int i=1; i<=j-1; ++i) {
                        h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
                        H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
                        H(i,j) =    h_ij;
                    }
                    nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
                    if (nu!=0.0) {
                        V[j] = (1. / H(j+1,j)) * w_j;
                    }
                    if (nu!=0.0) {
                        s(j) = H(j+1,j)/nu;
                        c(j) = H(j,j)/nu;
                        H(j,j) = nu;
                        H(j+1,j) = 0.0;
                        g(j+1) = -s(j)*g(j);
                        g(j) = c(j)*g(j);
                    }
                    rho = fabs(g(j+1));
                    std::cerr << "      iteration " << j << ", rho = " << rho << ", (tol = " << tol << ")" << std::endl;
                }

                if (j>=1) {
                    DeVector y(j);
                    for (int i=j; i>=1; --i) {
                        y(i) = g(i) / H(i,i);
                        for (int l=j; l>i; --l) {
                            y(i) -= H(i,l) * y(l) / H(i,i);
                        }
                    }
                    //std::cerr << "y = " << y << std::endl;
                    for (int i=1; i<=j; ++i) {
                        u += V[i-1] * y(i);
                    }
                    if (rho <= tol) {
                        std::cerr << "      iteration " << j << " return " << rho << std::endl;
                        res = rho;
                        return (iter-1)*CGMRES_ITER + j + 1;
                    }
                }
            }
            return maxIterations;
        }
        else {
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
            if (assemble_matrix==2) {
                A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
            }
            else {
                toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
            }

            if (Lambda.size() > 0) {
                DenseVector<Array<T> > rhs(N), x(N), residual(N), Ax(N);
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
                residual = Ax-rhs;
                res = std::sqrt(residual*residual);
                row_count = 1;
                for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
                    u[*row] = x(row_count);
                }
                return number_of_iterations;
            }
            else return -1;
        }
}

template <typename T, typename Index, typename MA>
int
GMRESM_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
             const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations, 
						 int assemble_matrix=1, int m=20)
{
       typedef typename IndexSet<Index >::const_iterator const_set_it;
       typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
       typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;
       typedef DenseVector<Array<T> >              DeVector;
       typedef GeMatrix<FullStorage<T, cxxblas::ColMajor> > DeMatrix;

       std::cerr << "GMRESM_Solve called." << std::endl;

       int N = Lambda.size();
       if (assemble_matrix==0) {
					std::cerr << " Algorithm GMRESM not implemented yet for assemble_matrix == 0" << std::endl;
					exit(1);
       }
       else {
           flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
           if (assemble_matrix==2) {
               A.toFlensSparseMatrix(Lambda, Lambda, A_flens, tol);
           }
           else {
               toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
           }

           if (Lambda.size() > 0) {
               DenseVector<Array<T> > rhs(N), x(N), residual(N), Ax(N);
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
               residual = Ax-rhs;
               res = std::sqrt(residual*residual);
               row_count = 1;
               for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
                   u[*row] = x(row_count);
               }
               return number_of_iterations;
           }
           else return -1;
       }							
}
						

template <typename T, typename Index, typename MA>
int
GMRES_Solve_PG(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol, MA &A, 
						Coefficients<Lexicographical,T,Index > &u, 
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
						int assemble_matrix = 1)
{
      typedef typename IndexSet<Index >::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
      typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;
      typedef DenseVector<Array<T> >              DeVector;
      typedef GeMatrix<FullStorage<T, cxxblas::ColMajor> > DeMatrix;

      std::cerr << "GMRES_Solve_PG called." << std::endl;

      int NumOfRows = (int)LambdaRow.size();
      int NumOfCols = (int)LambdaCol.size();
      if (assemble_matrix==0) {   //Algorithm by Saad, p.159
				std::cerr << " Algorithm GMRES_PG not implemented yet for assemble_matrix == 0" << std::endl;
				exit(1);      }
      else {
        flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
          if (assemble_matrix==2) {
						A.toFlensSparseMatrix(LambdaRow, LambdaCol, A_flens, tol);
          }
          else {
		        toFlensSparseMatrix(A, LambdaRow, LambdaCol, A_flens);
          }

          if (LambdaRow.size() > 0) {
              DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), residual(NumOfRows), Ax(NumOfRows);
              int row_count=1;
              const_coeff_it f_end = f.end();
              const_coeff_it u_end = u.end();
              for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
                  const_coeff_it f_it = f.find(*row);
                  if (f_it != f_end) rhs(row_count) = (*f_it).second;
                  else               rhs(row_count) = 0.;
                  const_coeff_it u_it = u.find(*row);
                  if (u_it != u_end) x(row_count) = (*u_it).second;
                  else               x(row_count) = 0.;

              }
              int number_of_iterations = lawa::gmres(A_flens,x,rhs, tol, maxIterations);
              Ax = A_flens*x;
              residual = Ax-rhs;
              res = std::sqrt(residual*residual);
              row_count = 1;
              for (const_set_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row, ++row_count) {
                  u[*row] = x(row_count);
              }
              return number_of_iterations;
          }
          else return -1;
      }							
							
}
						
						
template <typename T, typename Index, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRow, MA &A, Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
           int assemble_matrix)
{
        typedef typename IndexSet<Index >::const_iterator const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

        if (assemble_matrix==0) {   //Algorithm by Saad, p.237 ("CGNR")
            T alpha_cgls, beta_cgls, gammaPrev_cgls, gamma_cgls;
            Coefficients<Lexicographical,T,Index> r(SIZEHASHINDEX2D), q(SIZEHASHINDEX2D),
                                                  s(SIZEHASHINDEX2D), p(SIZEHASHINDEX2D);

            A.apply(u,0.,r,NoTrans);
            r -= f;
            r *= -1.;
            A.apply(r,0.,LambdaRow,s,Trans);
            p = s;
            gammaPrev_cgls = s*s;
            std::cerr << "      gammaPrev = " << gammaPrev_cgls << std::endl;

            for (int k=1; k<=maxIterations; k++) {
               q.setToZero();
               A.apply(p,0.,q,NoTrans);   //q = A*p;
               alpha_cgls = gammaPrev_cgls/(q*q);
               u +=   alpha_cgls*p;
               r -=   alpha_cgls*q;
               s.setToZero();
               A.apply(r,0.,LambdaRow,s,Trans);  // flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);

               gamma_cgls = s*s;
               res = r.norm(2.);
               std::cerr << "     iteration: " << k << " : current error ||A^T A u - A^f|| =" << sqrt(gamma_cgls)
                         << ", ||Au-f|| = " << res  << " (tol = " << tol << ")" << std::endl;
               if (sqrt(gamma_cgls)<=tol) {
                   Coefficients<Lexicographical,T,Index> help, AtAu;
                   A.apply(u, 0., help, cxxblas::NoTrans);
                   A.apply(help, 0., AtAu, cxxblas::Trans);
                   std::cerr << "      inner iterations: " << k  << ", ||A^T A u|| = " << AtAu.norm(2.) << std::endl;
                   return k-1;
               }

               beta_cgls  = gamma_cgls/gammaPrev_cgls;
               p *= beta_cgls;
               p += s;
               gammaPrev_cgls = gamma_cgls;
            }
        }
        else {
            //Attention: LambdaCol = supp(u)!!
            IndexSet<Index> LambdaCol;
            LambdaCol = supp(u);
            int NumOfRows = (int)LambdaRow.size();
            int NumOfCols = (int)LambdaCol.size();
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
            if (assemble_matrix==2) {
                A.toFlensSparseMatrix(LambdaRow, LambdaCol, A_flens, tol);
            }
            else {
                toFlensSparseMatrix(A, LambdaRow, LambdaCol, A_flens);
            }

            if (LambdaRow.size() > 0) {
                DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), residual(NumOfRows), Ax(NumOfRows);
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
                residual= Ax-rhs;
                res= std::sqrt(residual*residual);
                int col_count = 1;
                for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
                    u[*col] = x(col_count);
                }
                return number_of_iterations;
            }
            else return -1;
        }

}


template <typename T, typename Index, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol,  MA &A,
 					 Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, 
					 T &res, T tol, int maxIterations, int assemble_matrix)
{
		if (assemble_matrix==0) {
			std::cout << "CGLS called with LambdaRow and LambdaCol not set up to use with hashmap structure" << std::endl;
			exit(1);
		}
		else{
      typedef typename IndexSet<Index >::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
      typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

      int NumOfRows = (int)LambdaRow.size();
      int NumOfCols = (int)LambdaCol.size();
      flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(NumOfRows,NumOfCols);
      if (assemble_matrix==2) {
          A.toFlensSparseMatrix(LambdaRow, LambdaCol, A_flens, tol);
      }
      else {
          toFlensSparseMatrix(A, LambdaRow, LambdaCol, A_flens);
      }
      if (LambdaRow.size() > 0) {
          DenseVector<Array<T> > rhs(NumOfRows), x(NumOfCols), residual(NumOfRows), Ax(NumOfRows);
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
          residual= Ax-rhs;
          res= std::sqrt(residual*residual);
          int col_count = 1;
          for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
              u[*col] = x(col_count);
          }
          return number_of_iterations;
      }
      else return -1;
		}
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

}   // namespace lawa
