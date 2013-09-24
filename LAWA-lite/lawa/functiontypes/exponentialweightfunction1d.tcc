namespace lawa {

template <typename T>
T
ExponentialWeightFunction1D<T>::eta = 0.;

template <typename T>
T
ExponentialWeightFunction1D<T>::x1 = 0.;

template <typename T>
T
ExponentialWeightFunction1D<T>::x2 = 0.;

template <typename T>
bool
ExponentialWeightFunction1D<T>::with_precPoints = false;

template <typename T>
flens::DenseVector<Array<T> >
ExponentialWeightFunction1D<T>::singularPoints(2);

template <typename T>
void
ExponentialWeightFunction1D<T>::setEta(T _eta)
{
    setSingularPoints();
    eta=_eta;
}

template <typename T>
void
ExponentialWeightFunction1D<T>::setSingularPoints()
{
    //std::cout << "sing_pts = " << sing_pts << std::endl;
    singularPoints = -1., 1.;
}

template <typename T>
void
ExponentialWeightFunction1D<T>::setPrecPoints(T _x1, T _x2)
{
    with_precPoints=true;
    x1=_x1;
    x2=_x2;
}

template <typename T>
T
ExponentialWeightFunction1D<T>::alpha(T x)
{
    if (fabs(x)>=1.) {
        return -2*eta*fabs(x);
    }
    else {
        return 2*eta*(1./16.)*(-5-15*x*x+5*std::pow(x,4.)-std::pow(x,6.));
    }
}

template <typename T>
T
ExponentialWeightFunction1D<T>::dalpha(T x)
{
    if (x>1.) {
        return -2*eta;
    }
    else if (x<-1.) {
        return 2*eta;
    }
    else {
        return 2*eta*(1./16.)*(-30*x+20*std::pow(x,3.)-6*std::pow(x,5.));
    }
}

template <typename T>
T
ExponentialWeightFunction1D<T>::weight(T x)
{
    return with_precPoints ? exp(alpha(x)-0.5*(alpha(x1)+alpha(x2))) : exp(alpha(x));
}

template <typename T>
T
ExponentialWeightFunction1D<T>::sqrtweight(T x)
{
    return with_precPoints ? exp(0.5*(alpha(x)-0.5*(alpha(x1)+alpha(x2)))) : exp(0.5*alpha(x));
}

template <typename T>
T
ExponentialWeightFunction1D<T>::dweight(T x)
{
    return with_precPoints ? dalpha(x)*exp(alpha(x)-0.5*(alpha(x1)+alpha(x2))) : dalpha(x)*exp(alpha(x));
}


}   // namespace lawa
