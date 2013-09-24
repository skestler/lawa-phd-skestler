namespace lawa {

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS,ThetaTimeStepLocalOperator>::ThetaTimeStepSeparableRHS
(Function<T> &_fct_t, SpatialRHS &_F_x, ThetaTimeStepLocalOperator &_thetaTimeStepLocalOperator)
: fct_t(_fct_t), F_x(_F_x), thetaTimeStepLocalOperator(_thetaTimeStepLocalOperator),
  discrete_timepoint(0.), theta(1.), timestep(0.1)
{

}

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
void
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS,ThetaTimeStepLocalOperator>::setThetaTimeStepParameters
(T _theta,T _timestep, T _discrete_timepoint, const Coefficients<Lexicographical,T,Index> &_u_k)
{
    theta = _theta;
    timestep = _timestep;
    discrete_timepoint = _discrete_timepoint;
    u_k.clear();
    size_t hms;
    hms = _u_k.bucket_count();
    u_k.Rehash(hms);
    for (const_coeff_it it=_u_k.begin(); it!=_u_k.end(); ++it) {
        u_k[(*it).first] = (*it).second;
    }
    propagated_u_k.clear();
}

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
T
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS,ThetaTimeStepLocalOperator>::operator()(T t, const Index &index)
{
    return fct_t(t) * F_x(index);
}

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
void
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS,ThetaTimeStepLocalOperator>::initializePropagation
(const Coefficients<Lexicographical,T,Index> &f)
{
    bool propagationIsAlreadyComputed = true;
    if (propagated_u_k.size()==0) {
        size_t hms;
        hms = f.bucket_count();
        propagated_u_k.Rehash(hms);
        propagationIsAlreadyComputed = false;
    }

    for (const_coeff_it it=f.begin(); it!=f.end(); ++it) {
        if (propagated_u_k.find((*it).first)==propagated_u_k.end()) {
            propagationIsAlreadyComputed = false;
            break;
        }
    }

    if (!propagationIsAlreadyComputed) {
        //std::cerr << "Computing propagation for f with size " << f.size() << std::endl;
        for (const_coeff_it it=f.begin(); it!=f.end(); ++it) {
            propagated_u_k[(*it).first] = 0.;
        }
        if (theta!=1.) {
            thetaTimeStepLocalOperator.evalA(u_k, propagated_u_k, "residual");
            propagated_u_k *= (theta-1.)*timestep;
        }
        thetaTimeStepLocalOperator.evalM(u_k, propagated_u_k, "residual");
    }
}

template <typename T, typename Index, typename SpatialRHS, typename ThetaTimeStepLocalOperator>
T
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS,ThetaTimeStepLocalOperator>::operator()(const Index &index)
{
    T spatial_val = F_x(index);
    coeff_it it = (propagated_u_k).find(index);

    T return_val =  timestep*(   fct_t(discrete_timepoint)         *   theta  *spatial_val );
    if (theta != 1.)  {
        return_val  += timestep*( fct_t(discrete_timepoint-timestep)*(1.-theta)*spatial_val );
    }
    if (it!=(propagated_u_k).end()) {
        return_val += (*it).second;
    }
    else {
        std::cerr << "Error: Value in propagated_u_k not found!!" << std::endl;
    }
    return return_val;
}

}   // namespace lawa
