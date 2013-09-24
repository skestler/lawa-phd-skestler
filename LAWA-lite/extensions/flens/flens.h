/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#include <extensions/flens/cg.h>
#include <extensions/flens/crs.h>
#include <extensions/flens/diagonalmatrix.h>
#include <extensions/flens/densefunctions.h>
#include <extensions/flens/gmres.h>
#include <extensions/flens/inv.h>
#include <extensions/flens/lapack.h>
#include <extensions/flens/lapack_flens.h>
#include <extensions/flens/lapack_generic.h>
#include <extensions/flens/sparsematrix.h>
#include <extensions/flens/sparse_blas.h>
#include <extensions/flens/sparse_blas_flens.h>
#include <extensions/flens/spy.h>
#include <extensions/flens/powermethod.h>
#include <extensions/flens/vectornorm.h>

