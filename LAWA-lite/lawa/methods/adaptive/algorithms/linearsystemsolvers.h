/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_LINEARSYSTEMSOLVERS_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_LINEARSYSTEMSOLVERS_H 1

#include <extensions/flens/cg.h>
#include <extensions/flens/gmres.h>
#include <lawa/flensforlawa.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>


namespace lawa {


/* @assemble_matrix: "0" -> operate only on hashmap data,
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
         const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
         T &timeMatrixVector, int assemble_matrix=1);

/* @assemble_matrix: "0" -> operate only on hashmap data,
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
            int assemble_matrix=1);

/* @assemble_matrix: "0" -> operate only on hashmap data,
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
GMRESM_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
             const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations, 
						 int assemble_matrix=1, int m=20);
						
/* @assemble_matrix: "0" -> operate only on hashmap data,
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
GMRES_Solve_PG(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, 
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
						int assemble_matrix = 1);

/* @assemble_matrix: "0" -> operate only on hashmap data,
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRow, MA &A, Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations,
           int assemble_matrix=1);

/* @assemble_matrix: "0" -> operate only on hashmap data, use version without LambdaCol for that!
                     "1" -> use standard routine to assemble stiffness matrix
                     "2" -> use operator intern routines to assemble stiffness matrix
 */
template <typename T, typename Index, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol,  MA &A,
 					 Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, 
					 T &res, T tol, int maxIterations, int assemble_matrix=1);

//todo: adapt and revise these methods!!!
template <typename T, typename Index, typename SpaceIndex, typename MA>
int
CGLS_Solve(const IndexSet<Index> &LambdaRowOp,
           const IndexSet<SpaceIndex> &LambdaRowInitCond, MA &A,
           const IndexSet<Index> &LambdaCol,
           Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f,
           const Coefficients<Lexicographical,T,SpaceIndex > &u0,
           T &res, T tol = 1e-6, int maxIterations = 1000);

}   // namespace lawa

#include <lawa/methods/adaptive/algorithms/linearsystemsolvers.tcc>

#endif  //  LAWA_METHODS_ADAPTIVE_ALGORITHMS_LINEARSYSTEMSOLVERS_H
