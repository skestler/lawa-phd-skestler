namespace lawa {


//--- function * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrate_f1(const IntegralExpWeight<Quad,First,Second> &integral)
{
    typedef typename First::T T;
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);

    Support<T> common = first.support(integral.j,integral.k);
    common.l1 *= integral.RightmLeft;
    common.l1 += integral.left;
    common.l2 *= integral.RightmLeft;
    common.l2 += integral.left;

    int p = integral.expweight.singularPoints.length();
    DenseVector<Array<T> > singularPoints;
    if (IsPrimal<First>::value or IsOrthogonal<First>::value) {
        if (p>0) {
            DenseVector<Array<T> > firstSingularPoints = first.singularSupport(integral.j1,integral.k1);

            firstSingularPoints *= integral.RightmLeft;
            firstSingularPoints += integral.left;

            int m = firstSingularPoints.length();
            singularPoints.engine().resize(m+p);

            std::merge(firstSingularPoints.engine().data(),
                       firstSingularPoints.engine().data() + m,
                       integral.expweight.singularPoints.engine().data(),
                       integral.expweight.singularPoints.engine().data() + p,
                       singularPoints.engine().data());
        } else {
            singularPoints = first.singularSupport(integral.j1,integral.k1);
        }
    } else {
        assert(0);
        return 0.;
    }

    typename First::T x1 = (first.support(integral.j1,integral.k1).l2+first.support(integral.j1,integral.k1).l1)/2.;
    typename First::T x2 = 0.;
    integral.expweight.setPrecPoints(x1,x2);

    T ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        T a = singularPoints(i);
        T b = singularPoints(i+1);
        if (a==b)             continue;
        if (b<=common.l1)      continue;
        else if (a>=common.l2) break;
        ret += integral.quadrature(singularPoints(i),singularPoints(i+1));
    }
    return ret;
}

//--- function * primal/dual/orthogonal * primal/dual/orthogonal)
template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrate_f2(const IntegralExpWeight<Quad,First,Second> &integral)
{
    typedef typename First::T T;
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);

    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
                 second.support(integral.j2,integral.k2),common)<=0) {
        return 0.;
    }
    common.l1 *= integral.RightmLeft;
    common.l1 += integral.left;
    common.l2 *= integral.RightmLeft;
    common.l2 += integral.left;

    DenseVector<Array<T> > singularPoints;
    int p = integral.expweight.singularPoints.length();

    if (!BothDual<First,Second>::value) { // we do have (additional) singular points
        // merge singular points of bsplines/wavelets and function to one list.
        // -> implicit assumption: singular points are sorted!
        if (BothPrimal<First,Second>::value or BothOrthogonal<First,Second>::value) {
            DenseVector<Array<T> > firstSingularPoints = first.singularSupport(integral.j1,integral.k1);
            DenseVector<Array<T> > secondSingularPoints = second.singularSupport(integral.j2,integral.k2);

            firstSingularPoints *= integral.RightmLeft;
            firstSingularPoints += integral.left;
            secondSingularPoints *= integral.RightmLeft;
            secondSingularPoints += integral.left;

            int m = firstSingularPoints.length();
            int n = secondSingularPoints.length();
            singularPoints.engine().resize(m+n+p);

            std::merge(firstSingularPoints.engine().data(),
                       firstSingularPoints.engine().data() + m,
                       secondSingularPoints.engine().data(),
                       secondSingularPoints.engine().data() + n,
                       singularPoints.engine().data());
            if (p>0) {
                singularPoints(_(m+n+1,m+n+p)) = integral.expweight.singularPoints;
                std::inplace_merge(singularPoints.engine().data(),
                                   singularPoints.engine().data() + m + n,
                                   singularPoints.engine().data() + m + n + p);
            }
        }
    }
    else {
        assert(0);
        return (T)0.;
    }

    typename First::T x1 = (first.support(integral.j1,integral.k1).l2+first.support(integral.j1,integral.k1).l1)/2.;
    typename Second::T x2 = (second.support(integral.j2,integral.k2).l2+second.support(integral.j2,integral.k2).l1)/2.;
    integral.expweight.setPrecPoints(x1,x2);

    //T ret = 0.0;
    long double ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        T a = singularPoints(i);
        T b = singularPoints(i+1);
        if (a==b)             continue;
        if (b<=common.l1)      continue;
        else if (a>=common.l2) break;
        else                   ret += integral.quadrature(a,b);
    }

    return ret;
}


//--- function * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDualOrOrthogonal<First>::value, typename First::T>::Type
_integrand_f1(const IntegralExpWeight<Quad,First,Second> &integral, typename First::T x)
{
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);

    return    integral.expweight.weight(x)
            * first((x-integral.left)/(integral.RightmLeft),
                    integral.j1,integral.k1,integral.deriv1)
            / (integral.SqrtRightmLeft*std::pow(integral.RightmLeft,integral.deriv1));
}

//--- function * primal/dual/orthogonal * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDualOrOrthogonal<First>::value
                    and PrimalOrDualOrOrthogonal<Second>::value, typename First::T>::Type
_integrand_f2(const IntegralExpWeight<Quad,First,Second> &integral, typename First::T x)
{
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);

    return    integral.expweight.weight(x)
            * first((x-integral.left)/(integral.RightmLeft),
                    integral.j1,integral.k1,integral.deriv1)
            * second((x-integral.left)/(integral.RightmLeft),
                     integral.j2,integral.k2,integral.deriv2)
            / (integral.RightmLeft*std::pow(integral.RightmLeft,integral.deriv1+integral.deriv2));
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralExpWeight<Quad,First,Second>::integrand(typename First::T x) const
{
   return _f2 ? _integrand_f2(*this,x) : _integrand_f1(*this,x);
}

template <QuadratureType Quad, typename First, typename Second>
IntegralExpWeight<Quad,First,Second>::IntegralExpWeight(const First &_first, const T eta,
                                                        const T _left, const T _right)
    : expweight(), first(_first), second(_first), left(_left), right(_right),
      RightmLeft(right-left), SqrtRightmLeft(std::sqrt(right-left)),
      OneDivSqrtRightmLeft(1./SqrtRightmLeft),quadrature(*this), _f2(false)
{
    expweight.setEta(eta);
}

template <QuadratureType Quad, typename First, typename Second>
IntegralExpWeight<Quad,First,Second>::IntegralExpWeight(const First &_first, const Second &_second,
                                                        const T eta, const T _left, const T _right)
    : expweight(), first(_first), second(_second), left(_left), right(_right),
      RightmLeft(right-left), SqrtRightmLeft(std::sqrt(right-left)),
      OneDivSqrtRightmLeft(1./SqrtRightmLeft), quadrature(*this), _f2(true)

{
    expweight.setEta(eta);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralExpWeight<Quad,First,Second>::operator()(int _j, long _k, XType _e, int _deriv) const
{
    j1 = _j; k1 = _k; e1 = _e; deriv1 = _deriv;
    return _integrate_f1(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralExpWeight<Quad,First,Second>::operator()(int _j1, long _k1, XType _e1, int _deriv1,
                                                 int _j2, long _k2, XType _e2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; e2 = _e2; deriv2 = _deriv2;
    return _integrate_f2(*this);
}


}   // namespace lawa
