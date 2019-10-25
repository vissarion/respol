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

#ifndef HASHED_DETERMINANT_CGAL_2_H
#define HASHED_DETERMINANT_CGAL_2_H

#include "hashed_determinant_base.h"

#include <CGAL/Cartesian_d.h>

template <class _NT>
class HashedDeterminantCGAL2:public HashedDeterminantBase<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminantBase<NT>               HD;
        typedef typename HD::Index                      Index;
        typedef CGAL::Cartesian_d<NT>                   CK;
        typedef typename CK::LA                         LA;
        public:
        HashedDeterminantCGAL2():HD(){};
        HashedDeterminantCGAL2(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantCGAL2(Iterator begin,Iterator end):HD(begin,end){};
        const char* algorithm()const{return "CGAL_2\0";};
        inline NT compute_determinant(const Index&);
};

template <class _NT>
inline
typename HashedDeterminantCGAL2<_NT>::NT
HashedDeterminantCGAL2<_NT>::compute_determinant(
                const typename HashedDeterminantCGAL2<_NT>::Index &idx){
        int d = idx.size();
        //std::cout << first-last << "|" << d << std::endl;
        typename LA::Matrix A(d/2);
        typename LA::Matrix B(d/2);
        typename LA::Matrix C(d/2);
        typename LA::Matrix D(d/2);
        //std::vector<CPoint_d>::iterator s = first;

        for( int j = d/2; j < d; ++j ){
                //std::cout << *s << std::endl;
                for( int i = 0; i < d/2; ++i ){
                        std::cout << i << "," << j;
                        A(i,j-d/2) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = d/2; j < d; ++j ){
                //std::cout << *s << std::endl;
                for( int i = d/2; i < d; ++i ){
                        std::cout << i << "," << j;
                        B(i-d/2,j-d/2) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = 0; j < d/2; ++j ){
                //std::cout << *s << std::endl;
                for( int i = 0; i < d/2; ++i ){
                        std::cout << i << "," << j;
                        C(i,j) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = 0; j < d/2; ++j ){
                //std::cout << *s << std::endl;
                for( int i = d/2; i < d; ++i ){
                        std::cout << i << "," << j;
                        A(i-d/2,j) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        exit(0);
        return LA::determinant(M);
}

#endif // HASHED_DETERMINANT_CGAL_2_H
