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

#ifndef HASHED_DETERMINANT_BIRD_H
#define HASHED_DETERMINANT_BIRD_H

#include "hashed_determinant.h"
#include "bird.h"

template <class _NT>
class HashedDeterminantBird:public HashedDeterminant<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminant<NT>                   HD;
        typedef typename HD::Index                      Index;
        public:
        HashedDeterminantBird():HD(){};
        HashedDeterminantBird(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantBird(Iterator begin,Iterator end):HD(begin,end){}
        const char* algorithm()const{return "Bird\0";};
        inline NT compute_determinant(const Index&);
};

template <class _NT>
inline
typename HashedDeterminantBird<_NT>::NT
HashedDeterminantBird<_NT>::compute_determinant(
                const typename HashedDeterminantBird<_NT>::Index &idx){
        int d=idx.size();
        svector<svector<NT> > M;
        M.reserve(d);
        svector<NT> thisrow;
        thisrow.reserve(d);
        for(int j=0;j<d;++j){
                thisrow.clear();
                for(int i=0;i<d;++i)
                        thisrow.push_back(this->_points[idx[i]][j]);
                M.push_back(thisrow);
        }
        return Bird::determinant(M);
}

#endif // HASHED_DETERMINANT_BIRD_H
