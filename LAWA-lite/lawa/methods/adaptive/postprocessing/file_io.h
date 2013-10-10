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

#ifndef LAWA_METHODS_ADAPTIVE_POSTPROCESSING_FILE_IO_H
#define LAWA_METHODS_ADAPTIVE_POSTPROCESSING_FILE_IO_H 1

#include <iostream>
#include <fstream>
#include <cstring>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template<typename T, typename Index>
void
writeCoefficientsToFile(Coefficients<Lexicographical,T,Index> &u, int i, const char* filename);


template<typename T, typename Index>
void
readCoefficientsFromFile(Coefficients<Lexicographical,T,Index2D> &u, const char* filename);


/*
 * Required file format (example):
 * scaling,0,0,scaling,0,0 0.10241655476041190698
 * scaling,0,0,scaling,0,0 0.5411640725468976898
 * ...
 * No newlines at the beginning, indices are separated by kommas, the value at the end by a blank.
 */
template<typename T>
void
readCoefficientsFromFile(Coefficients<Lexicographical,T,Index2D> &u, const char* filename);

/*
 * Required file format (example):
 * scaling,0,0,scaling,0,0,scaling,0,0 0.10241655476041190698
 * scaling,0,0,scaling,0,0,scaling,0,1 0.5411640725468976898
 * ...
 * No newlines at the beginning, indices are separated by kommas, the value at the end by a blank.
 */
template<typename T>
void
readCoefficientsFromFile(Coefficients<Lexicographical,T,Index3D> &u, const char* filename);

}   // namespace lawa

#include <lawa/methods/adaptive/postprocessing/file_io.tcc>

#endif // LAWA_METHODS_ADAPTIVE_POSTPROCESSING_FILE_IO_H
