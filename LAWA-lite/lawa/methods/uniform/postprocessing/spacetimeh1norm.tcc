#include <iostream>

namespace lawa {

template<typename T, typename Basis2D>
SpaceTimeH1Norm<T,Basis2D>::SpaceTimeH1Norm(const Basis2D& _basis)
    : basis(_basis), integral_t(_basis.first, _basis.first), integral_x(_basis.second, _basis.second)
{
}

template<typename T, typename Basis2D>
T
SpaceTimeH1Norm<T,Basis2D>::operator()(XType xtype_t, int j_t, int k_t, 
                                       XType xtype_x, int j_x, int k_x) const
{
    T norm_t_L2 = integral_t(j_t, k_t, xtype_t, 0, j_t, k_t, xtype_t, 0);
    T norm_x_H1 = integral_x(j_x, k_x, xtype_x, 0, j_x, k_x, xtype_x, 0) + integral_x(j_x, k_x, xtype_x, 1, j_x, k_x, xtype_x, 1);

    // norm = norm_t_L2 * norm_x_H1 + norm_t_H1 * norm_x_H-1
    // Use:
    //      (norm_t_H1)^2 ~ 2^(2j) (norm_t_L2)^2
    //      (norm_x_H1')^2 ~ (norm_x_H1)^(-2)      
    return std::sqrt(norm_t_L2 * norm_x_H1 + lawa::pow2i<T>(2*j_t) * norm_t_L2 / norm_x_H1);                              
}

} // namespace lawa

