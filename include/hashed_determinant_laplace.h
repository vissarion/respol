// Copyright 2011-2012 National and Kapodistrian University of Athens,
// Greece.
//
// This file is part of respol.
//
// Respol is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// Respol is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with respol, see
// <http://www.gnu.org/licenses/>.

#ifndef HASHED_DETERMINANT_LAPLACE_H
#define HASHED_DETERMINANT_LAPLACE_H

#include "hashed_determinant_base.h"

template <class _NT>
class HashedDeterminantLaplace:public HashedDeterminantBase<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminantBase<NT>               HD;
        public:
        HashedDeterminantLaplace():HD(){};
        HashedDeterminantLaplace(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantLaplace(Iterator begin,Iterator end):HD(begin,end){}
        const char* algorithm()const{return "Laplace\0";};
};

#endif // HASHED_DETERMINANT_LAPLACE_H
