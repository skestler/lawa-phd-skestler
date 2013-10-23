namespace lawa {

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec)
: basis(_basis), Op(_Op), F(_F), Prec(_Prec), /*f_eps(_f_eps),*/ IsMW(false), sparsetree(false),
  alpha(0.5), gamma(0.1), residualType("standard"), hashMapSize(SIZEHASHINDEX2D),
  compute_f_minus_Au_error(false), write_coefficients_to_file(false)
{

}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::setParameters
(T _alpha, T _gamma, const char* _residualType, const char* _treeType, bool _IsMW,
 /*bool _compute_f_minus_Au_error,*/ bool _write_coefficients_to_file, size_t _hashMapSize)
{
    alpha = _alpha;
    gamma = _gamma;
    residualType = _residualType;
    if (strcmp(_treeType,"sparsetree")==0) sparsetree = true;
    else                                   sparsetree = false;
    IsMW = _IsMW;
    hashMapSize = _hashMapSize;
//    compute_f_minus_Au_error = _compute_f_minus_Au_error;
    write_coefficients_to_file = _write_coefficients_to_file;
    std::cerr << "Residual Type: " << residualType << std::endl;
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
cg_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, int NumOfIterations, T _init_cgtol,
         T EnergyNormSquared, const char *filename, const char *coefffilename, int maxDof)
{
    Coefficients<Lexicographical,T,Index> r(hashMapSize),       // residual vector for cg
                                          p(hashMapSize),       // auxiliary vector for cg
                                          Ap(hashMapSize);      // auxiliary vector for cg
    Coefficients<Lexicographical,T,Index> res(hashMapSize);     // approximate residual for f-Au
    Coefficients<Lexicographical,T,Index> u_leafs(hashMapSize); // "leafs" of u

    Index maxIndex;
    Index maxWaveletIndex;
    int *jmax = new int[1];
    int arrayLength = 1;
    long double Residual = 1.L;

    Timer time;
    Timer iteration_time;
    Timer galerkin_time;

    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        u_leafs[(*it).first] = 0.;
        r[(*it).first] = 0.;
        p[(*it).first] = 0.;
        Ap[(*it).first] = 0.;
        res[(*it).first] = 0.;
    }

    std::ofstream file(filename);
    file.precision(16);

    T time_total_comp = 0.;
    T time_total_galerkin = 0.;
    T time_total_residual = 0.;
    for (int iter=1; iter<=NumOfIterations; ++iter) {
        T time_mv_linsys = 0.;
        T time_mv_residual = 0.;
        T time_multitree_residual = 0.;
        iteration_time.start();

        int N = u.size();
        int N_residual = 0;


        std::cerr << std::endl << "   *****  AWGM Iteration " << iter << " *****" << std::endl << std::endl;
        std::cerr << "      Current number of dof = " << u.size() << std::endl;

        /* ******************* Resetting of vectors *********************** */

        //std::cerr << "      Resetting vectors..." << std::endl;

        //F.initializePropagation(u);
        std::cerr << "Attention: Propagation for finance problems not activated." << std::endl;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            T tmp = Prec((*it).first) * F((*it).first);
            r[(*it).first] = 0.;
            p[(*it).first] = tmp;
            Ap[(*it).first] = 0.;
        }

        std::cerr.precision(16);
        std::cerr << "            || f ||_{ell_2} = " << p.norm((T)2.) << std::endl;
        std::cerr.precision(6);
        //std::cerr << "      ...finished." << std::endl;

        /* ******************* CG method for Galerkin system *********************** */


        T time_galerkin = 0.;
        galerkin_time.start();
        //std::cerr << "      CG method started." << std::endl;
        //std::cerr << "DEBUG: Size of r = " << r.size() << std::endl;
        Op.eval(u,r,Prec,"galerkin");
        //std::cerr << "DEBUG: Size of r = " << r.size() << std::endl;
        r -= p;
        p = r;
        p *= (T)(-1.);
        T cg_rNormSquare = r*r;
        T tol = std::min(_init_cgtol,gamma*(T)Residual);
        int maxIterations=1000;
        int cg_iter=0;
        for (cg_iter=0; cg_iter<maxIterations; ++cg_iter) {
            if (std::sqrt(cg_rNormSquare)<=tol) {
                std::cerr << "         CG stopped with error " << sqrt(cg_rNormSquare) << std::endl;
                break;
            }
            //std::cerr << "    Iteration " << cg_iter+1 << std::endl;
            Ap.setToZero();
            time.start();
            Op.eval(p,Ap,Prec,"galerkin");
            time.stop();
            //std::cerr << "      DEBUG: " << time.elapsed() << std::endl;
            time_mv_linsys += time.elapsed();
            T pAp = p * Ap;
            T alpha = cg_rNormSquare/pAp;
            p *= alpha;
            u += p;
            p *= (T)1./alpha;
            Ap *= alpha;
            r += Ap;

            T cg_rNormSquarePrev = cg_rNormSquare;
            cg_rNormSquare = r*r;
            //std::cerr << "            Current error in cg: " << std::sqrt(cg_rNormSquare) << std::endl;
            //std::cerr << "            ||u||_2 = " << u.norm(2.) <<  std::endl;
            T beta = cg_rNormSquare/cg_rNormSquarePrev;
            p *= beta;
            p -= r;
        }
        std::cerr << "      CG method finished after " << cg_iter << " iterations." << std::endl;
        if (cg_iter>=maxIterations) {
            std::cerr << "      CG method exceeded maximum number of iterations " << maxIterations << std::endl;
        }
        time_mv_linsys *= 1./cg_iter;
        if (write_coefficients_to_file) {
            writeCoefficientsToFile(u, iter, coefffilename);
        }

        galerkin_time.stop();
        time_galerkin = galerkin_time.elapsed();
        time_total_galerkin += time_galerkin;
        /* ************************************************************************* */

        if (NumOfIterations==1)     return 0.;

        if (u.size()>=maxDof)       return Residual;

        /* ***************** Errors and some info on the solution ******************** */
        T EnergyError = 0., f_minus_Au_error = 0.;
        /*
        std::cerr << "      Computing energy error..." << std::endl;
        Ap.setToZero();
        A.eval(u,Ap,Prec,"galerkin");
        T uAu = Ap*u;
        T fu  = 0.;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            fu += (*it).second*Prec((*it).first) * F((*it).first);
        }
        EnergyError =  sqrt(fabs(EnergyNormSquared - 2*fu + uAu));
        std::cerr << "   ... finished with energy error: " << EnergyError << std::endl;

        if (compute_f_minus_Au_error) {
            compute_f_minus_Au(u, 1e-9, f_minus_Au_error);
        }
        */
        /*
        getLevelInfo(u, maxIndex, maxWaveletIndex, jmax, arrayLength);
        int J = -100;
        for (int i=0; i<arrayLength; ++i) {
            if (J<jmax[i]) J=jmax[i];
        }

        std::cerr << "   Level information: " << std::endl;
        std::cerr << "      maxIndex = " << maxIndex << ", maxWaveletIndex = "
                 << maxWaveletIndex << std::endl;
        std::cerr << "      arrayLength = " << arrayLength << std::endl;
        std::cerr << "      Highest level per coordinate direction:";
        for (int i=0; i<arrayLength; ++i) {
            std::cerr << " " << jmax[i];
        }
        std::cerr << std::endl;
        */


        //bool useSupportCenter=true;
        //std::stringstream coefffilename;
        //coefffilename << "u_coeff_mw_awgm_poisson3d_" << iter << ".dat";
        //plotScatterCoeff(u, basis, coefffilename.str().c_str(), useSupportCenter);

        //std::cerr << "Please hit enter." << std::endl;
        //getchar();


        /* ************************************************************************* */

        /* ******************* Computing approximate residual ********************** */

        //std::cerr << "      Computing residual..." << std::endl;
        //std::cerr << "         Computing multi-tree for residual..." << std::endl;
        //std::cerr << "           #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        time.start();
        res.setToZero();
        extendMultiTree(basis, u_leafs, res, residualType, IsMW, sparsetree);
        //extendMultiTreeAtBoundary(basis, u, res, J+1, sparsetree);
        time.stop();
        time_multitree_residual = time.elapsed();
        N_residual = res.size();
        //std::cerr << "         ... finished after " << time.elapsed() << std::endl;
        //std::cerr << "      #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        //std::cerr << "         Computing matrix vector product..." << std::endl;
        time.start();
        //A.eval(u,res,Prec);

        //F.initializePropagation(res);
        std::cerr << "Attention: Propagation for finance problems not activated." << std::endl;

        Op.eval(u,res,Prec,"residual");
        //Op.eval(u,res,Prec,"residual_standard");
        //std::cerr << "         ... finished." << std::endl;
        //std::cerr << "         Substracting right-hand side..." << std::endl;
        for (coeff_it it=res.begin(); it!=res.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        Residual = res.norm(2.);
        time_mv_residual = time.elapsed();
        time_multitree_residual += time_mv_residual;
        time_total_residual += time_multitree_residual;
        //std::cerr << "      ... finished." << std::endl;
        std::cerr << "      Residual: " << Residual << std::endl;
        if (Residual <= _eps) {
            std::cerr << "      Target tolerance reached: Residual = " << Residual
                      << ", eps = " << _eps << std::endl;
            return Residual;
        }

        /* ************************************************************************* */




        /* ********************** Computing next index set  ************************ */
        long double P_Lambda_Residual_square = 0.0L;
        if (u.size() > 0) {
            for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
                P_Lambda_Residual_square += std::pow(r[(*it).first],(T)2.);
                res.erase((*it).first);
            }
        }
        if (res.size()!=0) {
            T threshbound = std::sqrt(1-alpha*alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
            Coefficients<Bucket,T,Index> r_bucket;
            r_bucket.bucketsort(res, threshbound);
            //std::cerr << "         ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_Residual_square)
            //          << ", alpha*Residual = " << alpha*Residual << std::endl;

            for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
                P_Lambda_Residual_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
                r_bucket.addBucketToCoefficients(p,i);
                if (P_Lambda_Residual_square >= alpha*Residual*alpha*Residual) {
                    //r_bucket.addBucketToCoefficients(p,i+1);
                    break;
                }
            }
        }
        // Above we set res = res-res|_{supp u}. Now we change this back.
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            res[(*it).first] = 0.;
        }

        // The vector p satisfies the bulk criterion. But it is not yet a multi-tree...
        //std::cerr << "      Size of u before extension: " << u.size() << std::endl;
        //std::cerr << "      Size of extension set: " << p.size() << std::endl;
        for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
            if (u.find((*it).first)==u.end()) {
                completeMultiTree(basis, (*it).first, u, 0, sparsetree);
            }
        }
        //std::cerr << "      Size of u after extension: " << u.size() << std::endl;

        // Note that r is still supported on the previous Galerkin index set!
        u_leafs.clear();
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            if (r.find((*it).first)==r.end()) u_leafs[(*it).first] = 0.;
        }
        /* ************************************************************************* */

        iteration_time.stop();
        time_total_comp += iteration_time.elapsed();

        file << N << " " << N_residual << " " << time_total_comp << " " << EnergyError
                  << " " << Residual << " " << f_minus_Au_error << " "
                  << time_mv_linsys << " " << time_mv_residual << " " << time_multitree_residual
                  << " " << time_galerkin << " " << time_total_galerkin
                  << " " << time_total_residual << std::endl;



    }
    return Residual;
    file.close();
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
bicgstab_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, int NumOfIterations, T _init_cgtol,
               T EnergyNormSquared, const char *filename, const char *coefffilename, int maxDof)
{
    Coefficients<Lexicographical,T,Index> r(hashMapSize),       // residual vector for bicg
                                          rstar(hashMapSize),   // residual vector for bigcg
                                          p(hashMapSize),       // auxiliary vector for bicg
                                          s(hashMapSize),       // auxiliary vector for bicg
                                          Ap(hashMapSize),      // auxiliary vector for bicg
                                          As(hashMapSize);      // auxiliary vector for bicg
    Coefficients<Lexicographical,T,Index> res(hashMapSize);     // approximate residual for f-Au
    Coefficients<Lexicographical,T,Index> u_leafs(hashMapSize); // "leafs" of u

    Index maxIndex;
    Index maxWaveletIndex;
    int *jmax = new int[1];
    int arrayLength = 1;
    long double Residual = 1.L;

    Timer time;
    Timer iteration_time;
    Timer galerkin_time;

    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        u_leafs[(*it).first] = 0.;
        r[(*it).first] = 0.;
        rstar[(*it).first] = 0.;
        p[(*it).first] = 0.;
        s[(*it).first] = 0.;
        Ap[(*it).first] = 0.;
        As[(*it).first] = 0.;
        res[(*it).first] = 0.;
    }

    std::ofstream file(filename);
    file.precision(16);

    T time_total_comp = 0.;
    T time_total_galerkin = 0.;
    T time_total_residual = 0.;
    for (int iter=1; iter<=NumOfIterations; ++iter) {
        T time_mv_linsys = 0.;
        T time_mv_residual = 0.;
        T time_multitree_residual = 0.;
        iteration_time.start();

        int N = u.size();
        int N_residual = 0;


        std::cerr << std::endl << "   *****  AWGM Iteration " << iter << " *****" << std::endl << std::endl;
        std::cerr << "      Current number of dof = " << u.size() << std::endl;

        /* ******************* Resetting of vectors *********************** */

        //std::cerr << "      Resetting vectors..." << std::endl;

        F.initializePropagation(u);
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            T tmp = Prec((*it).first) * F((*it).first);
            r[(*it).first] = 0.;
            rstar[(*it).first] = 0.;
            p[(*it).first] = tmp;
            s[(*it).first] = 0.;
            Ap[(*it).first] = 0.;
            As[(*it).first] = 0.;
        }

        std::cerr.precision(16);
        std::cerr << "            || f ||_{ell_2} = " << p.norm((T)2.) << std::endl;
        std::cerr.precision(6);
        //std::cerr << "      ...finished." << std::endl;

        /* ******************* BiCGStab method for Galerkin system *********************** */


        T time_galerkin = 0.;
        galerkin_time.start();
        //std::cerr << "      CG method started." << std::endl;
        //std::cerr << "DEBUG: Size of r = " << r.size() << std::endl;
        Op.eval(u,r,Prec,"galerkin");
        //std::cerr << "DEBUG: Size of r = " << r.size() << std::endl;
        r     -= p;         // r = Ax - b
        r     *= (-1.);
        rstar  = r;
        p = r;

        T bicg_rNormSquare = r*r;
        T tol = std::min(_init_cgtol,gamma*(T)Residual);
        int maxIterations=1000;
        int bicg_iter=0;
        for (bicg_iter=0; bicg_iter<maxIterations; ++bicg_iter) {
            if (std::sqrt(bicg_rNormSquare)<=tol) {
                std::cerr << "         BiCGStab stopped with error " << sqrt(bicg_rNormSquare) << std::endl;
                break;
            }
            Ap.setToZero();
            time.start();
            Op.eval(p,Ap,Prec,"galerkin");
            time.stop();
            time_mv_linsys += time.elapsed();
            T r_times_rstar  = r*rstar;
            T rstar_times_Ap = rstar*Ap;
            T alpha = r_times_rstar / rstar_times_Ap;

            s  = Ap;
            s *= (-alpha);
            s += r;

            As.setToZero();
            Op.eval(s,As,Prec,"galerkin");
            T s_times_As  = s*As;
            T As_times_As = As*As;
            T omega = s_times_As / As_times_As;

            p *= alpha;
            s *= omega;
            u += p;
            u += s;
            p *= (T)(1./alpha);
            s *= (T)(1./omega);

            As *= (-omega);
            r   = As;
            r  += s;

            T r_times_rstar_next = r*rstar;
            T beta = (r_times_rstar_next / r_times_rstar) * (alpha / omega);
            Ap *= (-omega);
            p  += Ap;
            p  *= beta;
            p  += r;

            bicg_rNormSquare = r*r;
            std::cerr << "            Current error in bicg: " << std::sqrt(bicg_rNormSquare) << std::endl;
            std::cerr << "            ||u||_2 = " << u.norm(2.) <<  std::endl;
        }
        std::cerr << "      BiCGStab method finished after " << bicg_iter << " iterations." << std::endl;
        if (bicg_iter>=maxIterations) {
            std::cerr << "      BiCGStab method exceeded maximum number of iterations " << maxIterations << std::endl;
        }
        time_mv_linsys *= 1./bicg_iter;
        if (write_coefficients_to_file) {
            writeCoefficientsToFile(u, iter, coefffilename);
        }

        galerkin_time.stop();
        time_galerkin = galerkin_time.elapsed();
        time_total_galerkin += time_galerkin;
        /* ************************************************************************* */

        if (NumOfIterations==iter)     return 0.;

        if (u.size()>=maxDof)       return Residual;

        /* ***************** Errors and some info on the solution ******************** */
        T EnergyError = 0., f_minus_Au_error = 0.;



        /* ************************************************************************* */

        /* ******************* Computing approximate residual ********************** */

        //std::cerr << "      Computing residual..." << std::endl;
        //std::cerr << "         Computing multi-tree for residual..." << std::endl;
        //std::cerr << "           #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        time.start();
        res.setToZero();
        extendMultiTree(basis, u_leafs, res, residualType, IsMW, sparsetree);
        //extendMultiTreeAtBoundary(basis, u, res, J+1, sparsetree);
        time.stop();
        time_multitree_residual = time.elapsed();
        N_residual = res.size();
        //std::cerr << "         ... finished after " << time.elapsed() << std::endl;
        //std::cerr << "      #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        //std::cerr << "         Computing matrix vector product..." << std::endl;
        time.start();
        //A.eval(u,res,Prec);
        F.initializePropagation(res);
        Op.eval(u,res,Prec,"residual");
        //Op.eval(u,res,Prec,"residual_standard");
        //std::cerr << "         ... finished." << std::endl;
        //std::cerr << "         Substracting right-hand side..." << std::endl;
        for (coeff_it it=res.begin(); it!=res.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        Residual = res.norm(2.);
        time_mv_residual = time.elapsed();
        time_multitree_residual += time_mv_residual;
        time_total_residual += time_multitree_residual;
        //std::cerr << "      ... finished." << std::endl;
        std::cerr << "      Residual: " << Residual << std::endl;
        if (Residual <= _eps) {
            std::cerr << "      Target tolerance reached: Residual = " << Residual
                      << ", eps = " << _eps << std::endl;
            return Residual;
        }

        /* ************************************************************************* */




        /* ********************** Computing next index set  ************************ */
        long double P_Lambda_Residual_square = 0.0L;
        if (u.size() > 0) {
            for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
                P_Lambda_Residual_square += std::pow(r[(*it).first],(T)2.);
                res.erase((*it).first);
            }
        }
        if (res.size()!=0) {
            T threshbound = std::sqrt(1-alpha*alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
            Coefficients<Bucket,T,Index> r_bucket;
            r_bucket.bucketsort(res, threshbound);
            //std::cerr << "         ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_Residual_square)
            //          << ", alpha*Residual = " << alpha*Residual << std::endl;

            for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
                P_Lambda_Residual_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
                r_bucket.addBucketToCoefficients(p,i);
                if (P_Lambda_Residual_square >= alpha*Residual*alpha*Residual) {
                    //r_bucket.addBucketToCoefficients(p,i+1);
                    break;
                }
            }
        }
        // Above we set res = res-res|_{supp u}. Now we change this back.
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            res[(*it).first] = 0.;
        }

        // The vector p satisfies the bulk criterion. But it is not yet a multi-tree...
        //std::cerr << "      Size of u before extension: " << u.size() << std::endl;
        //std::cerr << "      Size of extension set: " << p.size() << std::endl;
        for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
            if (u.find((*it).first)==u.end()) {
                completeMultiTree(basis, (*it).first, u, 0, sparsetree);
            }
        }
        //std::cerr << "      Size of u after extension: " << u.size() << std::endl;

        // Note that r is still supported on the previous Galerkin index set!
        u_leafs.clear();
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            if (r.find((*it).first)==r.end()) u_leafs[(*it).first] = 0.;
        }
        /* ************************************************************************* */

        iteration_time.stop();
        time_total_comp += iteration_time.elapsed();

        file << N << " " << N_residual << " " << time_total_comp << " " << EnergyError
                  << " " << Residual << " " << f_minus_Au_error << " "
                  << time_mv_linsys << " " << time_mv_residual << " " << time_multitree_residual
                  << " " << time_galerkin << " " << time_total_galerkin
                  << " " << time_total_residual << std::endl;



    }
    return Residual;
    file.close();
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
approxL2(Coefficients<Lexicographical,T,Index> &u, T _eps, T (*weightFunction)(const Index &index),
         int NumOfIterations, T _init_cgtol,
         T EnergyNormSquared, const char *filename, const char *coefffilename)
{
    if (!IsMW) {    // Local operator is assumed to be a mass matrix
        this->cg_solve(u, _eps, NumOfIterations, _init_cgtol, EnergyNormSquared, filename, coefffilename);
    }
    else {          // Mass matrix is identity
        Coefficients<Lexicographical,T,Index> res(hashMapSize);
        Coefficients<Lexicographical,T,Index> p(hashMapSize);
        Coefficients<Lexicographical,T,Index> u_leafs(hashMapSize); // "leafs" of u

        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            u[(*it).first] = F((*it).first);
            u_leafs[(*it).first] = 0.;
            p[(*it).first] = 0.;
            res[(*it).first] = 0.;
        }

        if (NumOfIterations==1) return;

        std::ofstream convfile("conv_u0.txt");

        for (int iter=1; iter<=NumOfIterations; ++iter) {
            int N = u.size();
            int N_residual = 0;


            std::cerr << std::endl << "   *****  AWGM L2-Approx Iteration " << iter << " *****" << std::endl << std::endl;
            std::cerr << "      Current number of dof = " << u.size() << std::endl;

            /* ******************* Resetting of vectors *********************** */

            res.setToZero();
            p.clear();
            extendMultiTree(basis, u_leafs, res, residualType, IsMW, sparsetree);
            std::cerr << "      #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
            std::cerr << "        Computing residual..." << std::endl;

            // Store current multi-tree structure of u
            u_leafs = u;
            //std::cerr << res << std::endl;
            for (coeff_it it=res.begin(); it!=res.end(); ++it) {
                if (u.find((*it).first)!=u.end()) {
                    res.erase((*it).first);
                }
                else {
                    //std::cerr << (*it).first << " : " << F((*it).first)  << " "
                    //          << weightFunction((*it).first) << " "
                    //          << F((*it).first) * weightFunction((*it).first) << std::endl;
                    (*it).second = F((*it).first) * weightFunction((*it).first);
                }
            }
            T Residual = res.norm(2.);
            std::cerr << "      ... finished with Residual: " << Residual << std::endl;
            if (Residual <= _eps) {
                std::cerr << "      Target tolerance reached: Residual = " << Residual
                          << ", eps = " << _eps << std::endl;
                return;
            }
            convfile << u.size() << " " << Residual << std::endl;

            /* ************************************************************************* */

            /* ********************** Computing next index set  ************************ */
            long double P_Lambda_Residual_square = 0.0L;
            if (res.size()!=0) {
                T threshbound = std::sqrt(1-alpha*alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
                Coefficients<Bucket,T,Index> r_bucket;
                r_bucket.bucketsort(res, threshbound);
                std::cerr << "         ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_Residual_square)
                          << ", alpha*Residual = " << alpha*Residual << std::endl;

                for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
                    P_Lambda_Residual_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
                    r_bucket.addBucketToCoefficients(p,i);
                    if (P_Lambda_Residual_square >= alpha*Residual*alpha*Residual) {
                        //r_bucket.addBucketToCoefficients(p,i+1);
                        break;
                    }
                }
            }

            // Above we set res = res-res|_{supp u}. Now we change this back.
            for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
                res[(*it).first] = (*it).second;
            }

            // The vector p satisfies the bulk criterion. But it is not yet a multi-tree...
            std::cerr << "      Size of u before extension: " << u.size() << std::endl;
            for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
                if (u.find((*it).first)==u.end()) {
                    completeMultiTree(basis, (*it).first, u, 0, sparsetree);
                }
            }
            std::cerr << "      Size of u after extension: " << u.size() << std::endl;

            // Note that u_leafs is still supported on the previous Galerkin index set!
            for (coeff_it it=u.begin(); it!=u.end(); ++it) {
                if (u_leafs.find((*it).first)!=u_leafs.end()) u_leafs.erase((*it).first);
                else {
                    (*it).second = F((*it).first);//res[(*it).first];
                    u_leafs[(*it).first] = 0.;
                }
            }
            /* ************************************************************************* */
        }
    }
}


/*
template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
compute_f_minus_Au(Coefficients<Lexicographical,T,Index> &u, T eps, T &f_minus_Au_error)
{
    Coefficients<Lexicographical,T,Index> Au(hashMapSize);
    std::cerr << "      APPLY started..." << std::endl;
    Op.apply(u,Au,Prec,eps);
    std::cerr << "      ... finished." << std::endl;
    Au -= f_eps;
    f_minus_Au_error = Au.norm(2.);
    //Au.clear();
    std::cerr << "DEBUG: #buckets in Ap = " << Au.bucket_count() << std::endl;
    return;// Au.norm(2.);
}
*/

}   // namespace lawa


/*
        AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Interval,Multi,Orthogonal,Interval,Multi> adHeOp2D(basis,0.);
        typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;

        time.start();
        IndexSet<Index2D> Lambda;
        Lambda = supp(u);
        SparseMatrixT SpaMatA(Lambda.size(),Lambda.size());
        adHeOp2D.toFlensSparseMatrix(Lambda, Lambda, SpaMatA, 5);
        time.stop();
        std::cerr << "   Required time for assembling: " << time.elapsed() << std::endl;
*/

/* ******************* Computing approximate residual ********************** */
   // Counter example residual evaluation:
   // P_{\tilde \Lambda^{(1)}_k} (A \otimes I \cdots I) w_{\Lambda_k} = P_{\tilde \Lambda_k} (A \otimes I \cdots I) w_{\Lambda_k}
   // does not hold!
        /*
        std::cerr << "   Computing residual..." << std::endl;
        res.setToZero();
        std::cerr << "     Computing multi-tree for residual..." << std::endl;
        std::cerr << "        #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        time.start();
        Coefficients<Lexicographical,T,Index3D> r_approx(hashMapSize), tmp;
        r_approx = u;
        r_approx.setToZero();
        extendMultiTree(basis, u, r_approx, residualType, IsMW);
        std::cerr << "      Size of r_approx: " << r_approx.size() << std::endl;
        A.eval(u,r_approx,Prec,1);
        tmp = r_approx;
        res += r_approx;

        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            if (    (*it).first.index2.xtype==XWavelet && (*it).first.index2.j==0 && (*it).first.index2.k==4
                 && (*it).first.index3.xtype==XBSpline && (*it).first.index3.j==0 && (*it).first.index3.k==5   ) {
                std::cerr << (*it).first << std::endl;
            }
            if (    (*it).first.index1.xtype==XWavelet && (*it).first.index1.j==1 && (*it).first.index1.k==11) {
                std::cerr << (*it).first << std::endl;
            }
        }

        r_approx.clear();
        r_approx = u;
        extendMultiTree(basis, u, r_approx, 1);
        A.eval(u,r_approx,Prec,1);
        tmp -= r_approx;
        for (const_coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
            if (fabs((*it).second)>1e-13) {
                std::cerr << "Difference: " << (*it).first << " " << (*it).second << std::endl;
                if (tmp.find((*it).first)!=tmp.end()) std::cerr << "Contained in reference res." << std::endl;
                if (r_approx.find((*it).first)!=r_approx.end()) std::cerr << "Contained in new res." << std::endl;
            }
        }


        r_approx.clear();
        r_approx = u;
        r_approx.setToZero();
        extendMultiTree(basis, u, r_approx, residualType, IsMW);
        std::cerr << "      Size of r_approx: " << r_approx.size() << std::endl;
        A.eval(u,r_approx,Prec,2);
        res += r_approx;

        r_approx.clear();
        r_approx = u;
        r_approx.setToZero();
        std::cerr << "      Size of r_approx: " << r_approx.size() << std::endl;
        extendMultiTree(basis, u, r_approx, residualType, IsMW);
        A.eval(u,r_approx,Prec,3);
        res += r_approx;
        N_residual = res.size();
        std::cerr << "     ... finished after " << time.elapsed() << std::endl;
        std::cerr << "   #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        std::cerr << "     Substracting right-hand side..." << std::endl;
        for (coeff_it it=res.begin(); it!=res.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        time_multitree_residual = time.elapsed();
        time_mv_residual = 0;
        Residual = res.norm(2.);
        std::cerr << "   ... finished with Residual: " << Residual << std::endl;
        */
        /* ************************************************************************* */
