namespace lawa {

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::SYM_APPLY_1D(const Parameters &_parameters,
                                                          const Basis1D &_basis, MA &_A)
: parameters(_parameters), basis(_basis), A(_A)
{
}

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index>
SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::operator()(const Coefficients<Lexicographical,T,Index> &v, int k)
{
    int d=basis.d;
    Coefficients<Lexicographical,T,Index> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index> temp;
    temp = v;

    if (parameters.w_XBSpline) {
        int s = 0, count = 0;
        for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index> Lambda_v;
            Lambda_v=lambdaTilde1d_PDE((*it).second, basis,(k-s), basis.j0, (*it).second.j+(k-s)+1,false);
            for (const_set_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
    }
    else {
        /*
        int s = 0, count = 0;
        for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index1D> Lambda_v;
            Lambda_v=lambdaTilde1d_PDE_WO_XBSpline((*it).second, A.a.basis, (k-s), -50, 30);
            for (const_set_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
        */

        int s = 0, count = 0;
        T beta1 = 1./(d-0.5); //gamma-1 im nenner kuerzt sich raus!!
        T beta2 = beta1 - 0.1;

        for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index> Lambda_v;
            int kj = k-s;
            int minlevel, maxlevel;
            int j=(*it).second.j;
            if (j>=0) {
                maxlevel = j+kj;
                minlevel = j-kj;
                if (minlevel<0) {
                    minlevel =  int(std::floor(beta2*j-beta1*kj));
                }
            }
            else {//j<0
                minlevel = int(std::floor(j-beta1*kj));
                maxlevel = int(std::ceil(j+beta1*kj));
                if (maxlevel>0) maxlevel = int(std::ceil((j+beta1*kj)/beta2));
            }

            Lambda_v=lambdaTilde1d_PDE_WO_XBSpline((*it).second, A.a.basis, kj, minlevel, maxlevel);
            for (const_set_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }

    }
    return ret;
}

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index>
SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::operator()(const Coefficients<Lexicographical,T,Index> &v, int k, int J)
{
    //int d=basis.d;
    Coefficients<Lexicographical,T,Index> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index> temp;
    temp = v;

    if (parameters.w_XBSpline) {
        int s = 0, count = 0;
        for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
            IndexSet<Index> Lambda_v;
            Lambda_v=lambdaTilde1d_PDE((*it).second, basis,(k-s), basis.j0, J, false);
            for (const_set_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*it).second) * (*it).first;
            }
            ++count;
            s = int(log(T(count))/log(T(2))) + 1;
        }
    }
    else {
        assert(0);

    }
    return ret;
}

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index>
SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::operator()(const Coefficients<Lexicographical,T,Index> &v, T eps) {
    Coefficients<AbsoluteValue,T,Index> v_abs;
    v_abs = v;
    int k = SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::findK(v_abs, eps);
    Coefficients<Lexicographical,T,Index> ret;
    ret = SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::operator()(v, k);
    return ret;
}

template <typename T, typename Index, typename Basis1D, typename Parameters, typename MA>
int
SYM_APPLY_1D<T,Index,Basis1D,Parameters,MA>::findK(const Coefficients<AbsoluteValue,T,Index> &v, T eps) {
    int d=basis.d;
    if (v.size() == 0) return 1;
    T s=d-1.5;    //s = gamma-1, gamma the smoothness index of the wavelet basis

    T tau = 1.0 / (s + 0.5);
    // here the constant in (7.27) (KU-Wavelet) is estimated with 10
    int k_eps = static_cast<int>(10*log(std::pow(eps, -1.0/s)*std::pow(v.wtauNorm(tau), 1.0/s)) / log(2.0));
    DenseVector<Array<T> > normsec = v.norm_sections();
    T ErrorEstimateFactor = 1.;
    //std::cout << "eps = " << eps << ", k_eps = " << k_eps << std::endl;

    for (int k=1; k<=k_eps; ++k) {
        //std::cout << "At k = " << setw(3) << k;

        T R_k = 0.0;
        for (int i=k; i<=normsec.lastIndex()-1; ++i) {
            R_k += normsec(i+1);
        }
        R_k *= parameters.CA;
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k += std::pow(2.,-k*s) * normsec(1);
        //std::cout << ", R_k = " << setw(10) << R_k;

        for (int l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()-1) {
                //R_k += std::pow(l,-1.01)*std::pow(2.,-l*s) * normsec(k-l+1);
                R_k += std::pow(2.,-l*s) * normsec(k-l+1);
            }
        }
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k *= ErrorEstimateFactor;
        //std::cout << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10) << eps << endl;

        if (R_k<=eps) {
            std::cout << "   findK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            int maxlevel=22;
            if (d==2)         {    maxlevel=25; }
            else if (d==3)    {   maxlevel=19; }    //for non-singular examples, also lower values are possible
            return std::min(k,maxlevel);
        }
    }
    return std::min(k_eps,25);    //higher level differences result in translation indices that cannot be stored in int.
}


}    //namespace lawa

