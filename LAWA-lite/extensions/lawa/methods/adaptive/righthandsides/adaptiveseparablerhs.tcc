namespace lawa {

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3, typename RHSIntegral_x4,typename RHSIntegral_x5>
AdaptiveSeparableRhs<T,Index,RHSIntegral_x1,RHSIntegral_x2,RHSIntegral_x3,RHSIntegral_x4,RHSIntegral_x5>::
AdaptiveSeparableRhs(const RHSIntegral_x1 &_integral_x1, Coefficients<Lexicographical,T,Index1D> &_f_data_x1,
                     const RHSIntegral_x2 &_integral_x2, Coefficients<Lexicographical,T,Index1D> &_f_data_x2)
: integral_x1(_integral_x1), f_data_x1(_f_data_x1),
  integral_x2(_integral_x2), f_data_x2(_f_data_x2),
  integral_x3(_integral_x2), f_data_x3(_f_data_x2),
  integral_x4(_integral_x2), f_data_x4(_f_data_x2),
  integral_x5(_integral_x2), f_data_x5(_f_data_x2)
{

}

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3, typename RHSIntegral_x4,typename RHSIntegral_x5>
AdaptiveSeparableRhs<T,Index,RHSIntegral_x1,RHSIntegral_x2,RHSIntegral_x3,RHSIntegral_x4,RHSIntegral_x5>::
AdaptiveSeparableRhs(const RHSIntegral_x1 &_integral_x1, Coefficients<Lexicographical,T,Index1D> &_f_data_x1,
                     const RHSIntegral_x2 &_integral_x2, Coefficients<Lexicographical,T,Index1D> &_f_data_x2,
                     const RHSIntegral_x3 &_integral_x3, Coefficients<Lexicographical,T,Index1D> &_f_data_x3)
: integral_x1(_integral_x1), f_data_x1(_f_data_x1),
  integral_x2(_integral_x2), f_data_x2(_f_data_x2),
  integral_x3(_integral_x3), f_data_x3(_f_data_x3),
  integral_x4(_integral_x2), f_data_x4(_f_data_x2),
  integral_x5(_integral_x2), f_data_x5(_f_data_x2)
{

}

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3, typename RHSIntegral_x4,typename RHSIntegral_x5>
T
AdaptiveSeparableRhs<T,Index,RHSIntegral_x1,RHSIntegral_x2,RHSIntegral_x3,RHSIntegral_x4,RHSIntegral_x5>::
operator()(const Index2D &index2d)
{
    T val_x1 = 0.;
    const_coeff1d_it it_x1=f_data_x1.find(index2d.index1);
    if (it_x1!=f_data_x1.end()) val_x1 = (*it_x1).second;
    else {
        val_x1 = integral_x1(index2d.index1);
        f_data_x1[index2d.index1] = val_x1;
    }
    T val_x2 = 0.;
    const_coeff1d_it it_x2=f_data_x2.find(index2d.index2);
    if (it_x2!=f_data_x2.end()) val_x2 = (*it_x2).second;
    else {
        val_x2 = integral_x2(index2d.index2);
        f_data_x2[index2d.index2] = val_x2;
    }
    return val_x1 * val_x2;
}

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3, typename RHSIntegral_x4,typename RHSIntegral_x5>
T
AdaptiveSeparableRhs<T,Index,RHSIntegral_x1,RHSIntegral_x2,RHSIntegral_x3,RHSIntegral_x4,RHSIntegral_x5>::
operator()(const Index3D &index3d)
{
    T val_x1 = 0.;
    const_coeff1d_it it_x1=f_data_x1.find(index3d.index1);
    if (it_x1!=f_data_x1.end()) val_x1 = (*it_x1).second;
    else {
        val_x1 = integral_x1(index3d.index1);
        f_data_x1[index3d.index1] = val_x1;
    }
    T val_x2 = 0.;
    const_coeff1d_it it_x2=f_data_x2.find(index3d.index2);
    if (it_x2!=f_data_x2.end()) val_x2 = (*it_x2).second;
    else {
        val_x2 = integral_x2(index3d.index2);
        f_data_x2[index3d.index2] = val_x2;
    }
    T val_x3 = 0.;
    const_coeff1d_it it_x3=f_data_x3.find(index3d.index3);
    if (it_x3!=f_data_x3.end()) val_x3 = (*it_x3).second;
    else {
        val_x3 = integral_x3(index3d.index3);
        f_data_x3[index3d.index3] = val_x3;
    }
    return val_x1 * val_x2 * val_x3;
}

template <typename T, typename Index, typename RHSIntegral_x1, typename RHSIntegral_x2,
          typename RHSIntegral_x3, typename RHSIntegral_x4,typename RHSIntegral_x5>
void
AdaptiveSeparableRhs<T,Index,RHSIntegral_x1,RHSIntegral_x2,RHSIntegral_x3,RHSIntegral_x4,RHSIntegral_x5>::
clear()
{
    f_data_x1.clear();
    f_data_x2.clear();
    f_data_x3.clear();
    f_data_x4.clear();
    f_data_x5.clear();
}

}   // namespace lawa
