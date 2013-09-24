/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

Index1D::Index1D(void)
: j(0), k(0), xtype(XBSpline)
{
}

Index1D::Index1D(int _j, long _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype)
{

}

Index1D::Index1D(const Index1D &index)
: j(index.j), k(index.k), xtype(index.xtype)
{
}

Index1D::~Index1D(void)
{
}

short
Index1D::levelSum(void) const
{
    return this->j;
}


std::ostream& operator<<(std::ostream &s, const Index1D &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling," << _i.j << "," << _i.k;
    } else {
        s << "wavelet," << _i.j << "," << _i.k;
    }
    return s;
}


Index2D::Index2D(void)
: index1(), index2()
{
}

Index2D::Index2D(const Index1D &_index1, const Index1D &_index2)
: index1(_index1), index2(_index2)
{
}

Index2D::~Index2D(void)
{
}

short
Index2D::levelSum(void) const
{
    return this->index1.j+this->index2.j;
}

std::ostream& operator<<(std::ostream &s, const Index2D &_i)
{
    s << _i.index1 << "," << _i.index2;
    return s;
}


Index3D::Index3D(void)
: index1(), index2(), index3()
{
}

Index3D::Index3D(const Index1D &_index1, const Index1D &_index2, const Index1D &_index3)
: index1(_index1), index2(_index2), index3(_index3)
{
}

Index3D::~Index3D(void)
{
}

short
Index3D::levelSum(void) const
{
    return this->index1.j+this->index2.j+this->index3.j;
}

std::ostream& operator<<(std::ostream &s, const Index3D &_i)
{
    s <<  _i.index1 << "," << _i.index2 << "," << _i.index3;
    return s;
}


} //namespace lawa
