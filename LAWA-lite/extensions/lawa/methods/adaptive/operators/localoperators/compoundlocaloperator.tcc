namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(FirstLocalOperator &_firstLocalOp, SecondLocalOperator &_secondLocalOp)
: numOfLocalOp(2),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_secondLocalOp), fourthLocalOp(_secondLocalOp)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(FirstLocalOperator &_firstLocalOp, SecondLocalOperator &_secondLocalOp,
                        ThirdLocalOperator &_thirdLocalOp)
: numOfLocalOp(3),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_thirdLocalOp), fourthLocalOp(_secondLocalOp)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(const Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av)
{
    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            break;
        case 3:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            thirdLocalOp.eval(v, Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
       const char* evalType)
{
    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av, evalType);
            secondLocalOp.eval(v, Av, evalType);
            break;
        case 3:
            firstLocalOp.eval( v, Av, evalType);
            secondLocalOp.eval(v, Av, evalType);
            thirdLocalOp.eval(v, Av, evalType);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
template <typename Preconditioner>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
       Preconditioner &P)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }

    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            break;
        case 3:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            thirdLocalOp.eval(v, Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
template <typename Preconditioner>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v,
       Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, int operatornumber)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }

    switch (operatornumber)
    {
        case 1:
            firstLocalOp.eval( v, Av);
            break;
        case 2:
            secondLocalOp.eval(v, Av);
            break;
        case 3:
            thirdLocalOp.eval(v, Av);
            break;
        case 4:
            fourthLocalOp.eval(v, Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator: Non-admissible operatornumber " << operatornumber
                      << std::endl;
            exit(1);
    }
    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
template <typename Preconditioner>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
       Preconditioner &P, const char* evalType)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }

    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av, evalType);
            secondLocalOp.eval(v, Av, evalType);
            break;
        case 3:
            firstLocalOp.eval( v, Av, evalType);
            secondLocalOp.eval(v, Av, evalType);
            thirdLocalOp.eval(v, Av, evalType);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
template <typename Preconditioner>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>::
apply(Coefficients<Lexicographical,T,Index> &v,
      Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, T eps)
{
    int d = firstLocalOp.testBasis_CoordX.d;
    if (d!=secondLocalOp.testBasis_CoordX.d) {
        std::cerr << "CompoundLocalOperator<...>::apply: Order of first local operator is not equal"
                  << " to order of second local operator: "
                  << firstLocalOp.testBasis_CoordX.d << " " << secondLocalOp.testBasis_CoordX.d
                  << std::endl;
    }

    //todo: CA should be a read-in parameter
    T CA = 2.;
    Coefficients<Bucket,T,Index> v_bucket;
    T tol = 0.5*eps/CA;     //works better than tol=0.5*eps and multiply delta by CA... why??
    v_bucket.bucketsort(v,tol);
    long double squared_v_norm = (long double)std::pow(v.norm(2.),(T)2.);
    long double squared_v_bucket_norm = 0.;
    T delta=0.;
    int l=0;
    int support_size_all_buckets=0;
    for (int i=0; i<(int)v_bucket.buckets.size(); ++i) {
        squared_v_bucket_norm += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
        T squared_delta = fabs(squared_v_norm - squared_v_bucket_norm);
        support_size_all_buckets += v_bucket.buckets[i].size();
        delta = std::sqrt(squared_delta);
        l = i+1;
        if (squared_delta<tol*tol) {
            break;
        }
    }

    std::cerr << "      DEBUG: eps = " << eps << ", tol = " << tol
              << ", delta = " << delta << ", length of sum_p w_p = " << support_size_all_buckets << std::endl;
    if (delta>tol || eps < 5e-7) {
        std::cerr << "Warning: delta no longer computable!" << std::endl;
        delta = eps/2.;
    }
    delta = tol;
    Timer time;
    T critTime = 0.;
    int critBucketsize = 0;
    int critBucketnum = 0;
    T gamma = 0.5;
    if (d>=3) gamma = 1.5;      // smoothness of multiwavelets minus half the order of the operator
                                // multiwavelets are only C^1!!

    for (int i=0; i<l; ++i) {
        time.start();
        Coefficients<Lexicographical,T,Index> w_p;
        v_bucket.addBucketToCoefficients(w_p,i);
        if (w_p.size()==0) continue;
        T numerator = w_p.norm(2.) * support_size_all_buckets;
        T denominator = w_p.size() * (eps-delta) / CA;
        int jp = (int)std::max(((std::log(numerator/denominator) / std::log(2.) ) / gamma )/*-1*/, (T)0.);
        //int jp = ceil(std::max((std::log(numerator/denominator) / std::log(2.) / gamma )/*-1*/, (T)0.));

        //if (w_p.size()>1000 && eps > 5e-10) {
        //    jp = std::max(0, jp-1);
        //}
        for (const_coeff_it it=w_p.begin(); it!=w_p.end(); ++it) {
            if (numOfLocalOp==2) {
                T prec_col_index = P[(*it).first];
                IndexSet<Index1D> Lambda;
                Index1D coordX_col_index;
                typename FirstLocalOperator::notCoordXIndex notcoordX_col_index;

                firstLocalOp.split((*it).first, coordX_col_index, notcoordX_col_index);
                int maxlevel1 = std::min(coordX_col_index.j+jp,25);
                Lambda=lambdaTilde1d_PDE(coordX_col_index, firstLocalOp.testBasis_CoordX, jp,
                                           firstLocalOp.trialBasis_CoordX.j0, maxlevel1,false);
                firstLocalOp.nonTreeEval(coordX_col_index, notcoordX_col_index,
                                         prec_col_index*(*it).second, Lambda, Av, eps);
                secondLocalOp.split((*it).first, coordX_col_index, notcoordX_col_index);
                int maxlevel2 = std::min(coordX_col_index.j+jp,25);
                Lambda=lambdaTilde1d_PDE(coordX_col_index, secondLocalOp.testBasis_CoordX,jp,
                                           secondLocalOp.testBasis_CoordX.j0, maxlevel2,false);
                secondLocalOp.nonTreeEval(coordX_col_index, notcoordX_col_index,
                                          prec_col_index*(*it).second, Lambda, Av, eps);
            }
            else if (numOfLocalOp==3) {
                T prec_col_index = P[(*it).first];
                IndexSet<Index1D> Lambda;
                Index1D coordX_col_index;
                typename FirstLocalOperator::notCoordXIndex notcoordX_col_index;

                firstLocalOp.split((*it).first, coordX_col_index, notcoordX_col_index);
                int maxlevel1 = std::min(coordX_col_index.j+jp,25);
                Lambda=lambdaTilde1d_PDE(coordX_col_index, firstLocalOp.testBasis_CoordX, jp,
                                           firstLocalOp.trialBasis_CoordX.j0, maxlevel1,false);
                firstLocalOp.nonTreeEval(coordX_col_index, notcoordX_col_index,
                                         prec_col_index*(*it).second, Lambda, Av, eps);

                secondLocalOp.split((*it).first, coordX_col_index, notcoordX_col_index);
                int maxlevel2 = std::min(coordX_col_index.j+jp,25);
                Lambda=lambdaTilde1d_PDE(coordX_col_index, secondLocalOp.testBasis_CoordX,jp,
                                           secondLocalOp.testBasis_CoordX.j0, maxlevel2,false);
                secondLocalOp.nonTreeEval(coordX_col_index, notcoordX_col_index,
                                          prec_col_index*(*it).second, Lambda, Av, eps);

                thirdLocalOp.split((*it).first, coordX_col_index, notcoordX_col_index);
                int maxlevel3 = std::min(coordX_col_index.j+jp,25);
                Lambda=lambdaTilde1d_PDE(coordX_col_index, thirdLocalOp.testBasis_CoordX,jp,
                                           thirdLocalOp.testBasis_CoordX.j0, maxlevel3,false);
                thirdLocalOp.nonTreeEval(coordX_col_index, notcoordX_col_index,
                                          prec_col_index*(*it).second, Lambda, Av, eps);

            }
            else {
                std::cerr << "CompoundLocalOperator<...>::apply not yet implement for n>3." << std::endl;
                exit(1);
            }
        }
        time.stop();
        if (time.elapsed()>critTime) {
            critTime = time.elapsed();
            critBucketsize = w_p.size();
            critBucketnum = i;
        }
        //std::cout << "Bucket " << i << ": #wp= " << w_p.size() << ", jp=" << jp
        //          << ", time: " << time.elapsed() << std::endl;
    }
    std::cerr << "   Critical Bucket: " << critBucketnum << ", size = " << critBucketsize << ", time = " << critTime << std::endl;
    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
}

}   // namespace lawa
