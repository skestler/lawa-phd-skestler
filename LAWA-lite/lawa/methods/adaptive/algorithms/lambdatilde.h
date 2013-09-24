#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_LAMBDATILDE_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_LAMBDATILDE_H 1

#include <lawa/settings/enum.h>
#include <lawa/constructions/constructions.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

// Realizations of lambdaTilde for different Basis
template <typename T, typename Basis>
    IndexSet<Index1D>
    lambdaTilde1d_PDE(const Index1D &lambda, const Basis &basis,
                      int s_tilde, int jmin, int jmax, bool update);

template <typename T>
    IndexSet<Index1D>
    lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Periodic,CDF> &basis,
                      int s_tilde, int jmin, int jmax, bool update);

template <typename T, Construction Cons>
    IndexSet<Index1D>
    lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis, 
                      int s_tilde, int jmin, int jmax, bool update);


template <typename T>
    IndexSet<Index1D>
    lambdaTilde1d_PDE_WO_XBSpline(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis,
                                  int s_tilde, int jmin, int jmax);
template <typename T>
    IndexSet<Index1D>
    lambdaTilde1d_PDE_WO_XBSpline(const Index1D &lambda, const Basis<T,Orthogonal,R,Multi> &basis,
                                  int s_tilde, int jmin, int jmax);
                                  
template <typename T, Construction Cons>
    IndexSet<Index1D>
    lambdaTilde1d_WeightedPDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis,
                              int s_tilde_level, int jmin, int jmax, int s_tilde_singsupp=-1);

template <typename T>
    IndexSet<Index1D>
    lambdaTilde1d_WeightedPDE(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis,
                              int s_tilde_level, int jmin, int jmax, int s_tilde_singsupp=-1);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/lambdatilde.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_LAMBDATILDE_H

