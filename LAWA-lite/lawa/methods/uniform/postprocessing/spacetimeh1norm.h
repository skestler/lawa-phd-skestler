#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H 1

namespace lawa {

template <typename T, typename Basis2D>
class SpaceTimeH1Norm
{
    private:
        typedef typename Basis2D::FirstBasisType Basis_t;
        typedef typename Basis2D::SecondBasisType Basis_x;

        Integral<Gauss, Basis_t, Basis_t> integral_t;
        Integral<Gauss, Basis_x, Basis_x> integral_x; 

    public:

        const Basis2D &basis;

        SpaceTimeH1Norm(const Basis2D &_basis);

        /* Calculates the H1(0,T; H1)-norm (or W(0,T)-norm) for u = u1(t)*u2(x),
         *  ||u|| =  sqrt( ||u1||_L2 * ||u2||_H1 + |u1|_H1 * ||u2||_(H1)' )
         */
        T
        operator()(XType xtype_t, int j_t, int k_t, 
                   XType xtype_x, int j_x, int k_x) const;
};

} // namespace lawa

#include <lawa/methods/uniform/postprocessing/spacetimeh1norm.tcc>

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H

