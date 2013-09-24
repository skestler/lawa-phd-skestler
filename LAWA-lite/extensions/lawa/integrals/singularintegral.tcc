namespace lawa {



template <typename SingularKernel, typename First, typename Second>
typename RestrictTo<!IsDual<First>::value and !IsDual<Second>::value, typename First::T>::Type
_integrate(const SingularIntegral<SingularKernel,First,Second> &singularintegral)
{


    typedef typename First::T T;
    // note: called phi_... but can also be a wavelet.
    Support<T> supp_x = singularintegral.first.generator(singularintegral.e1).support(singularintegral.j1,singularintegral.k1);
    supp_x.l1 *= singularintegral.RightmLeft; supp_x.l1 += singularintegral.left;
    supp_x.l2 *= singularintegral.RightmLeft; supp_x.l2 += singularintegral.left;

    Support<T> supp_y = singularintegral.second.generator(singularintegral.e2).support(singularintegral.j2,singularintegral.k2);
    supp_y.l1 *= singularintegral.RightmLeft; supp_y.l1 += singularintegral.left;
    supp_y.l2 *= singularintegral.RightmLeft; supp_y.l2 += singularintegral.left;

    //std::cerr << "Integrate: [" << supp_x.l1 << ", " <<supp_x.l2 << "], [" << supp_y.l1 << ", " << supp_y.l2 << "]" << std::endl;
    DenseVector<Array<T> > SingularPoints_x =
             singularintegral.first.generator(singularintegral.e1).singularSupport(singularintegral.j1,singularintegral.k1);
    SingularPoints_x *= singularintegral.RightmLeft;
    SingularPoints_x += singularintegral.left;

    DenseVector<Array<T> > SingularPoints_y =
            singularintegral.second.generator(singularintegral.e2).singularSupport(singularintegral.j2,singularintegral.k2);
    SingularPoints_y *= singularintegral.RightmLeft;
    SingularPoints_y += singularintegral.left;

    std::vector<T> AllSingularPoints_x, AllSingularPoints_y;
    for (int i=SingularPoints_x.firstIndex(); i<=SingularPoints_x.lastIndex(); ++i) {
        T x_i = SingularPoints_x(i);
        AllSingularPoints_x.push_back(x_i);
        if (supp_y.l1 <= x_i && x_i <= supp_y.l2) AllSingularPoints_y.push_back(x_i);
    }
    for (int i=SingularPoints_y.firstIndex(); i<=SingularPoints_y.lastIndex(); ++i) {
        T y_i = SingularPoints_y(i);
        AllSingularPoints_y.push_back(y_i);
        if (supp_x.l1 <= y_i && y_i <= supp_x.l2) AllSingularPoints_x.push_back(y_i);
    }

    sort(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
    typename std::vector<T>::iterator it_x;
    it_x = unique(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
    AllSingularPoints_x.resize(it_x - AllSingularPoints_x.begin());

    typename std::vector<T>::iterator it_y;
    sort(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
    it_y = unique(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
    AllSingularPoints_y.resize(it_y - AllSingularPoints_y.begin());


    typename std::vector<T>::const_iterator it_x1, it_x2, it_y1, it_y2;
    it_x1=AllSingularPoints_x.begin();
    it_x2=AllSingularPoints_x.begin(); ++it_x2;

    long double ret = 0.L;
    while (it_x2!=AllSingularPoints_x.end()) {
        long double a1 = (*it_x1), b1 = (*it_x2);
        it_y1=AllSingularPoints_y.begin();
        it_y2=AllSingularPoints_y.begin(); ++it_y2;
        while (it_y2!=AllSingularPoints_y.end()) {
            long double a2 = (*it_y1), b2 = (*it_y2);
            long double tmp = singularintegral.singularquadrature(a1, b1, a2, b2, (long double)1e-10);
            ret += tmp;
            //std::cerr << "   [" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]: " << tmp << std::endl;
            ++it_y1; ++it_y2;
        }
        ++it_x1; ++it_x2;
    }

    return (T)ret;
}

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
typename FirstPolynomial::T
_integrate(const SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial> &singularintegral,
           typename FirstPolynomial::T a1, typename FirstPolynomial::T b1,
           typename FirstPolynomial::T a2, typename FirstPolynomial::T b2)
{
    typedef typename FirstPolynomial::T T;

    if (a2>b1+1e-14 || a1>b2+1e-14) {
        //std::cerr << " No singularity" << std::endl;
        return singularintegral.singularquadrature(a1, b1, a2, b2, (long double)1e-10);
    }
    else if (fabs(a1-b2)<=1e-14 || fabs(b1-a2)<=1e-14) {
        //std::cerr << " Corner singularity" << std::endl;
        return singularintegral.singularquadrature(a1, b1, a2, b2, (long double)1e-10);
    }
    else {
        //std::cerr << " Diagonal singularity" << std::endl;
        if (fabs(a2-a1)<1e-14 && fabs(b2-b1)<1e-14) {  //quadratic domain with diagonal singularity
            //std::cerr << "   Quadratic domain" << std::endl;
            return singularintegral.singularquadrature(a1, b1, a2, b2, (long double)1e-10);
        }
        DenseVector<Array<T> > SingularPoints_x(2); SingularPoints_x = a1, b1;
        DenseVector<Array<T> > SingularPoints_y(2); SingularPoints_y = a2, b2;

        std::vector<T> AllSingularPoints_x, AllSingularPoints_y;
        for (int i=SingularPoints_x.firstIndex(); i<=SingularPoints_x.lastIndex(); ++i) {
            T x_i = SingularPoints_x(i);
            AllSingularPoints_x.push_back(x_i);
            if (a2 < x_i && x_i < b2) AllSingularPoints_y.push_back(x_i);
        }
        for (int i=SingularPoints_y.firstIndex(); i<=SingularPoints_y.lastIndex(); ++i) {
            T y_i = SingularPoints_y(i);
            AllSingularPoints_y.push_back(y_i);
            if (a1 <= y_i && y_i <= b1) AllSingularPoints_x.push_back(y_i);
        }

        sort(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
        typename std::vector<T>::iterator it_x;
        it_x = unique(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
        AllSingularPoints_x.resize(it_x - AllSingularPoints_x.begin());

        typename std::vector<T>::iterator it_y;
        sort(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
        it_y = unique(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
        AllSingularPoints_y.resize(it_y - AllSingularPoints_y.begin());

        typename std::vector<T>::const_iterator it_x1, it_x2, it_y1, it_y2;
        it_x1=AllSingularPoints_x.begin();
        it_x2=AllSingularPoints_x.begin(); ++it_x2;

        long double ret = 0.L;
        while (it_x2!=AllSingularPoints_x.end()) {
            long double _a1 = (*it_x1), _b1 = (*it_x2);
            it_y1=AllSingularPoints_y.begin();
            it_y2=AllSingularPoints_y.begin(); ++it_y2;
            while (it_y2!=AllSingularPoints_y.end()) {
                long double _a2 = (*it_y1), _b2 = (*it_y2);
                T tmp = singularintegral.singularquadrature(_a1, _b1, _a2, _b2, (long double)1e-10);
                ret += tmp;
                ++it_y1; ++it_y2;
            }
            ++it_x1; ++it_x2;
        }
        return (T)ret;
    }

}

template <typename SingularKernel, typename First, typename Second>
SingularIntegral<SingularKernel,First,Second>::SingularIntegral(const SingularKernel &_singularkernel,
                                                                const First &_first,
                                                                const Second &_second,
                                                                const T _left, const T _right)
    : singularkernel(_singularkernel), first(_first), second(_second), left(_left), right(_right),
      RightmLeft(right-left), SqrtRightmLeft(std::sqrt(right-left)),
      OneDivSqrtRightmLeft(1./SqrtRightmLeft),
      singularquadrature(*this)
{
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::p1(T x) const
{
    return first.generator(e1).operator()((x-left)/(RightmLeft),j1,k1,deriv1) / (SqrtRightmLeft*RightmLeft);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::p2(T x) const
{
    //return second.generator(e2).operator()(x,j2,k2,deriv2);
    return second.generator(e2).operator()((x-left)/(RightmLeft),j2,k2,deriv2) / (SqrtRightmLeft*RightmLeft);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::kernel(T x) const
{
    return singularkernel.operator()(x);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::operator()
                                               (int _j1, long _k1, XType _e1, int _deriv1,
                                                int _j2, long _k2, XType _e2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; e2 = _e2; deriv2 = _deriv2;
    return _integrate(*this);
}




template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial>::SingularIntegralPP
                                                           (const SingularKernel &_singularkernel,
                                                            const FirstPolynomial &_first,
                                                            const SecondPolynomial &_second)
    : singularkernel(_singularkernel), first(_first), second(_second), singularquadrature(*this)
{
}

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
typename FirstPolynomial::T
SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial>::p1(T x) const
{
   return first.operator()(x);
}

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
typename FirstPolynomial::T
SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial>::p2(T x) const
{
    return second.operator()(x);
}

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
typename FirstPolynomial::T
SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial>::kernel(T x) const
{
   return singularkernel.operator()(x);
}

template <typename SingularKernel, typename FirstPolynomial, typename SecondPolynomial>
typename FirstPolynomial::T
SingularIntegralPP<SingularKernel,FirstPolynomial,SecondPolynomial>::operator()(T a1, T b1, T a2, T b2) const

{
   return _integrate(*this, a1, b1, a2, b2);
}

}   // namespace lawa

