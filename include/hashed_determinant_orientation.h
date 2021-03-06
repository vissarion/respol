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

#ifndef HASHED_DETERMINANT_ORIENTATION_H
#define HASHED_DETERMINANT_ORIENTATION_H

#include "hashed_determinant_base.h"

#include <CGAL/Cartesian_d.h>

template <class _NT>
class HashedDeterminantOrientation:public HashedDeterminantBase<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminantBase<NT>               HD;
        typedef typename HD::Index                      Index;
        typedef CGAL::Cartesian_d<NT>                   CK;
        typedef typename CK::LA                         LA;
        public:
        HashedDeterminantOrientation():HD(){};
        HashedDeterminantOrientation(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantOrientation(Iterator begin,Iterator end):
                HD(begin,end){}
        const char* algorithm()const{return "Orientation_CGAL\0";};
        inline NT compute_determinant(const Index&);
};

template <class _NT>
inline
typename HashedDeterminantOrientation<_NT>::NT
HashedDeterminantOrientation<_NT>::compute_determinant(
                const typename HashedDeterminantOrientation<_NT>::Index &idx){
        int d=idx.size()-1;
        typename LA::Matrix M(d);
        for(int j=0;j<d-1;++j)
                for(int i=1;i<=d;++i)
                        M(i-1,j)=this->_points[idx[i]][j]-
                                 this->_points[idx[0]][j];
        return LA::determinant(M);
}

#endif // HASHED_DETERMINANT_ORIENTATION_H
