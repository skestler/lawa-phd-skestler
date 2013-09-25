namespace lawa {

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
PayoffIntegral2D<Quad,Basis,PayoffFunction>::PayoffIntegral2D
(const Basis &_basis, const PayoffFunction &_payofffunction,
 const T _left_x1, const T _right_x1, const T _left_x2, const T _right_x2,
 bool _useSpecialRefinement, T _maxRectangleLength, int _order)
    : basis(_basis), payofffunction(_payofffunction),
      highestOrderQuadrature(*this), highOrderQuadrature(*this), mediumOrderQuadrature(*this),
      lowOrderQuadrature(*this),
      left_x1(_left_x1), right_x1(_right_x1), left_x2(_left_x2), right_x2(_right_x2),
      RightmLeft_x1(right_x1-left_x1), SqrtRightmLeft_x1(std::sqrt(right_x1-left_x1)),
      RightmLeft_x2(right_x2-left_x2), SqrtRightmLeft_x2(std::sqrt(right_x2-left_x2)),
      useSpecialRefinement(_useSpecialRefinement), maxRectangleLength(_maxRectangleLength), order(_order)
{
    int minOrder = basis.first.d+2;
    highestOrderQuadrature.setOrder(order);
    highOrderQuadrature.setOrder(std::max(order-3,minOrder));
    mediumOrderQuadrature.setOrder(std::max(order-6,minOrder));
    lowOrderQuadrature.setOrder(std::max(order-9,minOrder));
}

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
typename Basis::T
PayoffIntegral2D<Quad,Basis,PayoffFunction>::integrand(T x1, T x2) const
{
    const typename Basis::FirstBasisType::BasisFunctionType  &first  = basis.first.generator(e1);
    const typename Basis::SecondBasisType::BasisFunctionType &second = basis.second.generator(e2);

    T val_x1 = (1./SqrtRightmLeft_x1) * first((x1-left_x1)/(RightmLeft_x1),j1,k1,0);
    T val_x2 = (1./SqrtRightmLeft_x2) * second((x2-left_x2)/(RightmLeft_x2),j2,k2,0);

    T ret = val_x1 * val_x2 * payofffunction.payoff(x1,x2);
    return ret;
}

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
typename Basis::T
PayoffIntegral2D<Quad,Basis,PayoffFunction>::operator()(const Index2D &index2d) const
{
    j1 = index2d.index1.j;
    j2 = index2d.index2.j;
    k1 = index2d.index1.k;
    k2 = index2d.index2.k;
    e1 = index2d.index1.xtype;
    e2 = index2d.index2.xtype;

    const typename Basis::FirstBasisType::BasisFunctionType  &first  = basis.first.generator(e1);
    const typename Basis::SecondBasisType::BasisFunctionType &second = basis.second.generator(e2);

    DenseVector<Array<T> > singularPoints_x1, singularPoints_x2;

    // Computing singular points in x1 direction
    DenseVectorT singsupp_x1 = first.singularSupport(j1,k1);
    singsupp_x1 *= RightmLeft_x1;
    singsupp_x1 += left_x1;

    int m_x1 = singsupp_x1.length();
    int p_x1 = payofffunction.singPts_x1.length();
    singularPoints_x1.engine().resize(m_x1+p_x1);

    std::merge(singsupp_x1.engine().data(),
               singsupp_x1.engine().data() + m_x1,
               payofffunction.singPts_x1.engine().data(),
               payofffunction.singPts_x1.engine().data() + p_x1,
               singularPoints_x1.engine().data());

    // Computing singular points in x2 direction
    DenseVectorT singsupp_x2 = second.singularSupport(j2,k2);
    singsupp_x2 *= RightmLeft_x2;
    singsupp_x2 += left_x2;

    int m_x2 = singsupp_x2.length();
    int p_x2 = payofffunction.singPts_x2.length();
    singularPoints_x2.engine().resize(m_x2+p_x2);

    std::merge(singsupp_x2.engine().data(),
               singsupp_x2.engine().data() + m_x2,
               payofffunction.singPts_x2.engine().data(),
               payofffunction.singPts_x2.engine().data() + p_x2,
               singularPoints_x2.engine().data());

    //std::cerr << "PayoffIntegral2D<Quad,Basis,PayoffFunction>::operator() index = " << index2d << std::endl;
    long double ret = 0.L;
    for (int i1=singularPoints_x1.firstIndex(); i1<singularPoints_x1.lastIndex(); ++i1) {
        T a1 = singularPoints_x1(i1), b1 = singularPoints_x1(i1+1);
        for (int i2=singularPoints_x2.firstIndex(); i2<singularPoints_x2.lastIndex(); ++i2) {
            T a2 = singularPoints_x2(i2), b2 = singularPoints_x2(i2+1);
            //std::cerr << "[" << a1 << "," << b1 << "],[" << a2 << "," << b2 << "]" << std::endl;
            ret += (long double) this->integrate(a1, b1, a2, b2);
        }
    }
    return (T)ret;
}

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
typename Basis::T
PayoffIntegral2D<Quad,Basis,PayoffFunction>::integrate(T a1, T b1, T a2, T b2) const
{
    //std::cerr << "Integrating ["  << a1 << "," << b1 << "],[" << a2 << "," << b2 << "]" << std::endl;
    // Integrating over [a1,b1] x [a2,b2]

    if (!useSpecialRefinement) {
        T h1 = (b1-a1);
        T h2 = (b2-a2);
        if (h1>maxRectangleLength) return integrate(a1, a1+h1/2., a2, b2) + integrate(a1+h1/2., b1, a2, b2);
        if (h2>maxRectangleLength) return integrate(a1, b1, a2, a2+h2/2.) + integrate(a1, b1, a2+h2/2., b2);

        return highestOrderQuadrature(a1, b1, a2, b2);
    }
    else {
        T val_a1a2 = payofffunction.payoff(a1,a2);
        T val_a1b2 = payofffunction.payoff(a1,b2);
        T val_b1a2 = payofffunction.payoff(b1,a2);
        T val_b1b2 = payofffunction.payoff(b1,b2);

        if (!payofffunction.isCritical(a1,b1,a2,b2) && val_a1a2==0 && val_a1b2==0 && val_b1a2==0 && val_b1b2==0) {
            // Payoff-function is constantly zero
            //std::cerr << "   => ZERO." << std::endl;
            //T tmp = _getOrderAndValue(a1, b1, a2, b2);
            //if (fabs(tmp)>0) std::cerr << "Here is something going wrong..." << std::endl;

            return 0.;
        }

        if (val_a1a2>0 && val_a1b2>0 && val_b1a2>0 && val_b1b2>0) {
            if ( fabs(b1-a1)<4*maxRectangleLength && fabs(b2-a2)<4*maxRectangleLength ) {
                //std::cerr << "   => KINK integrate" << std::endl;
                return _getOrderAndValue(a1, b1, a2, b2);
            }
            else {
                T h1 = (b1-a1);
                T h2 = (b2-a2);
                //std::cerr << "   => SUBDIVIDE" << std::endl;
                if (h1>h2) return    integrate(a1, a1+h1/2., a2, b2) + integrate(a1+h1/2., b1, a2, b2);
                else       return    integrate(a1, b1, a2, a2+h2/2.) + integrate(a1, b1, a2+h2/2., b2);
            }
        }
        else {
            // Integration over kink
            //return _getOrderAndValue(a1, b1, a2, b2);
            if ( fabs(b1-a1)<maxRectangleLength && fabs(b2-a2)<maxRectangleLength ) {
                //std::cerr << "   => KINK integrate" << std::endl;
                return _getOrderAndValue(a1, b1, a2, b2);
            }
            else {
                T h1 = (b1-a1);
                T h2 = (b2-a2);
                //std::cerr << "   => SUBDIVIDE" << std::endl;
                if (h1>h2) return    integrate(a1, a1+h1/2., a2, b2) + integrate(a1+h1/2., b1, a2, b2);
                else       return    integrate(a1, b1, a2, a2+h2/2.) + integrate(a1, b1, a2+h2/2., b2);
            }
        }
    }
}

template <QuadratureType Quad, typename Basis, typename PayoffFunction>
typename Basis::T
PayoffIntegral2D<Quad,Basis,PayoffFunction>::_getOrderAndValue(T a1, T b1, T a2, T b2) const
{
    if (Quad == FullGridGL) {
        T h1 = (b1-a1), h2 = (b2-a2);
        if (h1*h2<1e-4) {
            return highestOrderQuadrature(a1,b1,a2,b2);
        }
        else if (h1*h2<1e-3) {
            return highestOrderQuadrature(a1,b1,a2,b2);
        }
        else if (h1*h2<1e-2) {
            return highestOrderQuadrature(a1,b1,a2,b2);
        }
        else {
            return highestOrderQuadrature(a1,b1,a2,b2);
        }
    }
}

}   // namespace lawa
